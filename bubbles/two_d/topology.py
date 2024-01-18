"""
TODO: Problem mit der Reihenfolge: Eigentlich sollte das Hinzufügen der bubbles unabhängig
davon sein, allerdings gilt das nicht, wenn bubbles mal collisions erlauben oder nicht. 
Das ist trivial.

TODO: Topology2 / 3 does not work, the intersection points are not added correctly, 
check out intersection with domain holes!!!.
"""

from __future__ import annotations
from typing import Any, Callable
from pydantic import BaseModel
import shapely
from shapely.geometry import Polygon

from bubbles.two_d.hole import Rectangular

from shapely.affinity import translate as shapely_translate
from tqdm import tqdm


class Bubble(BaseModel):
    polygon: shapely.Polygon
    level: int
    is_hole: bool

    class Config:
        arbitrary_types_allowed = True

    def __hash__(self) -> int:
        return self.polygon.__hash__()


class Topology(object):
    """
    1. Add hole/inclusion q

    2. Find intersections I = p1,...,pn (elements from self.polygons)

    # update the shape of q constrained by higher level bubbles
    3. Find pj1, .., pjk from I such that all have higher level
    3.1 If q is completely contained within merge(pj1,...,pjk) -> ignore this hole
    3.2 modify q = q.difference(pj1, .., pjk) => q might be a multipolygon

    # update the lower level bubbles according to the (remaining) q
    4.  Find pj_1, .., pj_l from I such that all have lower level
        and modify pj_k = pj_k.difference(q) for all k = 1 .. l  => pj_k might be a multipolygon
    4.1 For each original pj_k that is fully contained in q: remove pj_k  ( maybe pj_k is empty set then ? )
    4.2 Else if pj_k is a Multipolygon then split it into multiple polygons and delete pj_k

    # At this point q might be a MultiPolygon
    5. Find pi1,...pim, subset of I, s.t. all have the same level as q
    5.1 Merge pi1,...pim and q. This is a (possible) multipolygon denoted as P.
    5.2 Delete pi1,...pim and add each subpolygon of P.
    5.3 If P is a multipolygon split into polygons.

    # TODO: Separate bubbles and domain logic. The domain could be interpreted as a bubble with level 0, and
    behave exactly the same way as the other bubbles. Meaning, it can grow or disappear. In a second
    step this behaviour could be regulated if not desired.

    # TODO: Make main add_bubble functions static.
    """

    DEFAULT_GRID_SIZE = 1e-15

    def __init__(self, polygon: shapely.Polygon, grid_size=None):
        self.bubbles = {Bubble(polygon=polygon, level=0, is_hole=False)}

        self.ref_domain = Bubble(polygon=polygon, level=0, is_hole=False)
        # TODO: Check if it is enough just to set the grid_size level for all new geometries,
        # and then leave the parameter out for the rest of the code.
        # E.g. doing: polygon = shapely.set_precision(polygon, grid_size)

        self.grid_size = grid_size if grid_size else Topology.DEFAULT_GRID_SIZE

        # level2is_hole is a dictionary that maps level to is_hole
        # Once a bubble is added to the topology, the level2is_hole dictionary is updated,
        # the first time a bubble with a new level is added, the is_hole value is set and
        # cannot be changed.
        self.__level2is_hole = {0: False}

    @property
    def level2is_hole(self) -> dict[int, bool]:
        """Return a dictionary that maps level to is_hole."""
        return self.__level2is_hole

    @property
    def levels(self) -> list[int]:
        """Return all currently present levels."""
        return sorted(set([p.level for p in self.bubbles]))

    def __add_bubble_merge(self, bubble: Bubble) -> bool:
        """
        Add a bubble to the topology.
        """
        polygon, level, is_hole = bubble.polygon, bubble.level, bubble.is_hole
        if not level > 0:
            raise ValueError("Level must be positive integer.")
        if level in self.level2is_hole and self.level2is_hole[level] != is_hole:
            raise ValueError("Bubbles with same level must have same is_hole value.")

        # Clip the bubble with respect to the reference domain
        polygon = polygon.intersection(
            self.ref_domain.polygon, grid_size=self.grid_size
        )
        if polygon.area == 0:
            # new_bubble is completely outside the reference domain and is ignored
            return False
        if isinstance(polygon, shapely.GeometryCollection):
            polygon = self.fix_geometryCollection(polygon)

        # Find intersecting polygons and partition by level relation
        lower_bubbles = set()
        higher_bubbles = set()
        equal_bubbles = set()
        for _bubble in self.bubbles:
            if _bubble.polygon.intersects(polygon):
                if level > _bubble.level:
                    lower_bubbles.add(_bubble)
                elif level < _bubble.level:
                    higher_bubbles.add(_bubble)
                else:
                    equal_bubbles.add(_bubble)
                    if not _bubble.is_hole == is_hole:
                        raise ValueError(
                            "Bubbles with same level must have same is_hole value."
                        )

        # Iterate over higher level bubbles and update new_bubble
        if len(higher_bubbles) > 0:
            higher_polygons_union = shapely.union_all(
                [b.polygon for b in higher_bubbles]
            )
            polygon = polygon.difference(
                higher_polygons_union, grid_size=self.grid_size
            )
            if not polygon.area > self.grid_size:
                # new_bubble is completely contained in higher level bubbles and is ignored
                return False
            # new_bubble might be a multipolygon or geometry collection
            if isinstance(polygon, shapely.GeometryCollection):
                polygon = self.fix_geometryCollection(polygon)

        # Iterate over lower level bubbles and update them
        for lower_bubble in lower_bubbles:
            p_diff = lower_bubble.polygon.difference(polygon, grid_size=self.grid_size)
            if p_diff.area > self.grid_size:
                if isinstance(p_diff, shapely.GeometryCollection):
                    p_diff = self.fix_geometryCollection(p_diff)

                if isinstance(p_diff, shapely.Polygon):
                    # Add intersection points to new_bubble
                    polygon = self.add_interserction_points(polygon, p_diff)
                    if isinstance(polygon, shapely.GeometryCollection):
                        polygon = self.fix_geometryCollection(polygon)

                    # Add modified lower level bubble
                    self.bubbles.add(
                        Bubble(
                            polygon=p_diff,
                            level=lower_bubble.level,
                            is_hole=lower_bubble.is_hole,
                        )
                    )
                else:
                    for p_sub in p_diff.geoms:
                        # Add intersection points to new_bubble
                        polygon = self.add_interserction_points(polygon, p_sub)
                        if isinstance(polygon, shapely.GeometryCollection):
                            polygon = self.fix_geometryCollection(polygon)

                        # Add modified lower level bubble
                        self.bubbles.add(
                            Bubble(
                                polygon=p_sub,
                                level=lower_bubble.level,
                                is_hole=lower_bubble.is_hole,
                            )
                        )

            # Remove original lower bubble
            self.bubbles.remove(lower_bubble)

        # Remove equal level intersecting bubbles as they are replaced by
        # their union with new_bubble.
        for bubble in equal_bubbles:
            self.bubbles.remove(bubble)

        polygon = shapely.union_all(
            [polygon] + [p.polygon for p in equal_bubbles], grid_size=self.grid_size
        )
        if isinstance(polygon, shapely.GeometryCollection):
            polygon = self.fix_geometryCollection(polygon)

        # Add polygon to self.bubbles
        if isinstance(polygon, shapely.Polygon):
            self.bubbles.add(Bubble(polygon=polygon, level=level, is_hole=is_hole))
        else:
            for q_sub in polygon.geoms:
                self.bubbles.add(Bubble(polygon=q_sub, level=level, is_hole=is_hole))

        self.__level2is_hole[level] = is_hole
        return True

    def __add_bubble_disjoint(
        self, bubble: Bubble, bubble_distance: float, domain_distance: float
    ) -> bool:
        """Add a bubble to the topology.

        The bubble is added only if it is disjoint to all other bubbles, having minimal distance `bubble_distance`.

        The bubble can intersect the domain boundary, if `domain_distance` is 0, otherwise it must be disjoint and have
        distance `domain_distance` to the domain boundary.

        """
        if not bubble.level > 0:
            raise ValueError("Level must be positive integer.")
        if (
            bubble.level in self.level2is_hole
            and self.level2is_hole[bubble.level] != bubble.is_hole
        ):
            raise ValueError("Bubbles with same level must have same is_hole value.")
        if bubble_distance < self.grid_size:
            raise ValueError("Minimal distance must be larger than grid_size.")
        if domain_distance > 0:
            if domain_distance < self.grid_size:
                raise ValueError("Minimal distance must be larger than grid_size.")

        polygon = bubble.polygon
        # Assert that the polygon intersects with the reference domain
        if not self.ref_domain.polygon.intersects(polygon):
            return False

        # Clip the bubble with respect to the reference domain OR assert
        # minimal distance to the reference domain boundary.
        if domain_distance == 0:
            polygon = bubble.polygon.intersection(
                self.ref_domain.polygon, grid_size=self.grid_size
            )
            if isinstance(polygon, shapely.GeometryCollection):
                polygon = self.fix_geometryCollection(polygon)
        else:
            if self.ref_domain.polygon.boundary.distance(polygon) < domain_distance:
                return False

        # Check distance to other bubbles
        for _bubble in self.bubbles:
            # Only iterate over non-boundary bubbles
            if _bubble.level > 0:
                if _bubble.polygon.distance(polygon) < bubble_distance:
                    return False

        # Add bubble
        self.bubbles.add(
            Bubble(
                polygon=polygon,
                level=bubble.level,
                is_hole=bubble.is_hole,
            )
        )
        self.__level2is_hole[bubble.level] = bubble.is_hole
        return True

    def add(self, bubble: Bubble, **kwargs) -> bool:
        """
        Add a bubble to the topology.
        TODO: Merge certain levels and others disjoint.
        kwargs: flag to decide if bubble is merged with others having the same level.
        list of levels to be "kollisions Frei" with,
        list of distances to each level in the list.
        """
        if "bubble_distance" in kwargs and "domain_distance" in kwargs:
            return self.__add_bubble_disjoint(bubble, **kwargs)
        return self.__add_bubble_merge(bubble)

    def fix_geometryCollection(
        self, geo: shapely.GeometryCollection
    ) -> shapely.Polygon | shapely.MultiPolygon:
        """
        Fix a GeometryCollection by merging all polygons and removing all non-polygons.
        """
        polygons = [p for p in geo.geoms if isinstance(p, shapely.Polygon)]
        if len(polygons) == 0:
            raise ValueError("geo does not contain any polygons.")
        elif len(polygons) == 1:
            return polygons[0]
        else:
            union = shapely.union_all(polygons)
            if not isinstance(union, (shapely.Polygon, shapely.MultiPolygon)):
                raise ValueError(
                    f"union is of type {type(union)} but should be Polygon or MultiPolygon."
                )
            else:
                return union

    def add_interserction_points(
        self, multi_poly: shapely.Polygon | shapely.MultiPolygon, poly: shapely.Polygon
    ) -> shapely.Polygon | shapely.MultiPolygon:
        """
        Add intersection points of poly with multi_poly to multi_poly.
        """
        union = poly.union(multi_poly, grid_size=self.grid_size)
        if not isinstance(union, (shapely.Polygon, shapely.MultiPolygon)):
            if isinstance(union, shapely.GeometryCollection):
                union = self.fix_geometryCollection(union)
            else:
                raise ValueError(f"union is of type {type(union)}")
        polygon = multi_poly.intersection(
            union,
            grid_size=self.grid_size,
        )
        if not isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon)):
            if isinstance(polygon, shapely.GeometryCollection):
                polygon = self.fix_geometryCollection(polygon)
            else:
                raise ValueError(f"polygon is of type {type(polygon)}")
        return polygon

    def plot(self) -> Any:
        """
        Plot the reference domain and the bubbles.
        """
        plot(self.bubbles, self.level2is_hole, self.ref_domain)


def plot(bubbles, level2is_hole, ref_domain) -> Any:
    """
    Plot the reference domain and the bubbles.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    bubbles = list(bubbles)
    bubbles.sort(key=lambda x: x.level)

    hole_boundary_color = "orange"
    ref_boundary_color = "black"
    ref_domain_color = "silver"
    ref_domain_hole_color = "white"

    # level2is_hole = self.level2is_hole TODO: Remove this

    # Filter out levels that are not represented by any bubble.
    # This is necessary as previous added bubbles might have been removed.
    present_levels = set([bubble.level for bubble in bubbles])
    level2is_hole = {
        lvl: is_hole for lvl, is_hole in level2is_hole.items() if lvl in present_levels
    }
    levels = set(level2is_hole.keys())

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = ref_domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            lvl2cl[lvl] = plt.cm.cool(norm(lvl))
    # Make a legend.
    handles = [
        Line2D([0], [0], label="Reference boundary", color=ref_boundary_color),
        Line2D([0], [0], label="Hole boundary", color=hole_boundary_color),
    ]
    for lvl, cl in lvl2cl.items():
        if lvl == 0:
            label_text = f"Domain lvl = {lvl}"
        else:
            label_text = (
                f"Hole lvl = {lvl}" if level2is_hole[lvl] else f"Inclusion lvl = {lvl}"
            )
        handles.append(mpatches.Patch(color=cl, label=label_text))
    plt.legend(handles=handles)

    # Plot the reference domain.
    x, y = ref_domain.polygon.exterior.xy
    plt.plot(x, y, color=ref_boundary_color)
    plt.fill(x, y, color=ref_domain_color)
    for pp in ref_domain.polygon.interiors:
        x, y = pp.xy
        plt.plot(x, y, color=ref_boundary_color)
        plt.fill(x, y, color=ref_domain_hole_color)

    # Plot the bubbles
    for bubble in bubbles:
        # Filter out reference domain
        if bubble.level != 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Re-Plot Bubbles without holes
    for bubble in bubbles:
        if len(bubble.polygon.interiors) == 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Show the plot
    plt.show()


def plot2(bubbles, level2is_hole) -> Any:
    """
    Plot the reference domain and the bubbles.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    bubbles.sort(key=lambda x: x.level)

    hole_boundary_color = "orange"
    ref_boundary_color = "black"
    ref_domain_color = "silver"
    ref_domain_hole_color = "white"

    # level2is_hole = self.level2is_hole TODO: Remove this

    # Filter out levels that are not represented by any bubble.
    # This is necessary as previous added bubbles might have been removed.
    present_levels = set([bubble.level for bubble in bubbles])
    level2is_hole = {
        lvl: is_hole for lvl, is_hole in level2is_hole.items() if lvl in present_levels
    }
    levels = set(level2is_hole.keys())

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = ref_domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            lvl2cl[lvl] = plt.cm.cool(norm(lvl))

    # Plot the bubbles
    for bubble in bubbles:
        x, y = bubble.polygon.exterior.xy
        plt.fill(x, y, color=lvl2cl[bubble.level])
        if len(bubble.polygon.interiors) > 0:
            for pp in bubble.polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color="white")

    # Re-Plot Bubbles without holes
    for bubble in bubbles:
        if len(bubble.polygon.interiors) == 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Show the plot
    plt.show()


class PeriodicTopology(Topology):
    """
    # TODO: Add periodic flag to add function.
    # (x,y, or both)
    """

    def __init__(self, midpoint, width, height, periodic_x, periodic_y, grid_size=None):
        polygon = shapely.Polygon(
            Rectangular(midpoint=midpoint, width=width, height=height).discretize_hole(
                refs=4
            )
        )
        super().__init__(polygon, grid_size)
        self.bounds = (
            midpoint[0] - width / 2,
            midpoint[1] - height / 2,
            midpoint[0] + width / 2,
            midpoint[1] + height / 2,
        )
        self.width = width
        self.height = height
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y

    def add(self, bubble: Bubble, **kwargs) -> bool:
        poly = self.clip(bubble.polygon)

        if poly.area == 0:
            return False
        if isinstance(poly, shapely.Polygon):
            super().add(
                Bubble(polygon=poly, level=bubble.level, is_hole=bubble.is_hole),
                **kwargs,
            )
        else:
            for sub_poly in poly.geoms:
                super().add(
                    bubble=Bubble(
                        polygon=sub_poly, level=bubble.level, is_hole=bubble.is_hole
                    ),
                    **kwargs,
                )

    def clip(self, polygon):
        """Clip polygon along periodic boundaries until it is inside the reference domain."""
        poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
        self_x_min, self_y_min, self_x_max, self_y_max = self.bounds

        if not isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon)):
            raise ValueError(f"polygon is of type {type(polygon)}")

        if self.ref_domain.polygon.intersects(polygon):
            if poly_x_min < self_x_min and self.periodic_x:
                _in, _out = polygon.intersection(
                    self.ref_domain.polygon
                ), polygon.difference(self.ref_domain.polygon)
                polygon = shapely.union_all(
                    [_in, shapely_translate(_out, xoff=self.width)]
                )
            elif poly_y_min < self_y_min and self.periodic_y:
                _in, _out = polygon.intersection(
                    self.ref_domain.polygon
                ), polygon.difference(self.ref_domain.polygon)
                polygon = shapely.union_all(
                    [_in, shapely_translate(_out, yoff=self.height)]
                )
            if poly_x_max > self_x_max and self.periodic_x:
                _in, _out = polygon.intersection(
                    self.ref_domain.polygon
                ), polygon.difference(self.ref_domain.polygon)
                polygon = shapely.union_all(
                    [_in, shapely_translate(_out, xoff=-self.width)]
                )
            elif poly_y_max > self_y_max and self.periodic_y:
                _in, _out = polygon.intersection(
                    self.ref_domain.polygon
                ), polygon.difference(self.ref_domain.polygon)
                polygon = shapely.union_all(
                    [_in, shapely_translate(_out, yoff=-self.height)]
                )
            else:
                return polygon
        else:
            return shapely.Polygon()
        return self.clip(polygon)


class Bubble_v2(BaseModel):
    polygon: shapely.Polygon | shapely.MultiPolygon
    level: int
    is_hole: bool

    class Config:
        arbitrary_types_allowed = True

    def __hash__(self) -> int:
        return self.polygon.__hash__()


class Topology_v2:
    # TODO: bubbles (maybe as dict level to multipolygon), All bubbles with the same level
    # are merged into one multipolygon.

    # TODO: idea, define cutout levels, when topology is passed to gmsh api.
    """
    1. Sample n different Bubbles (different levels), each bubble is a polygon or Multipolygon.
    2. write a dict level to multipolygon, where all bubbles with the same level are merged into one multipolygon.
    3. Clipping: Clip each multipolygon with respect to the reference domain.
    4. For each level (lower to higher), difference the higher Multipolygons.
    """

    def __init__(self, clip_polygon=None) -> None:
        self.bubbles = []
        self.__level2is_hole = dict()
        self.clip_polygon = clip_polygon
        self.grid_size = 1e-15

    @property
    def level2is_hole(self) -> dict[int, bool]:
        """Return a dictionary that maps level to is_hole."""
        return self.__level2is_hole

    def add(
        self,
        polygon: shapely.Polygon | shapely.MultiPolygon,
        level: int,
        is_hole: bool,
        clip: bool = False,
    ) -> bool:
        # TODO: Check if clip_level idea is good or not. That is clipping a bubble with respect to a certain level.
        # The clipping polygon changes over time.
        if clip:
            polygon = Topology_v2.intersection(
                polygon, self.clip_polygon, grid_size=self.grid_size
            )

        # Make sure the bubble is not too small
        if polygon.area < self.grid_size:
            return False

        new_bubbles, is_added = self.__add_bubble_merge(
            self.bubbles, bubble, grid_size=self.grid_size
        )
        if is_added:
            # Make sure that all multipolygons are split into polygons
            bubbles_refactored = []
            for bub in new_bubbles:
                if isinstance(bub.polygon, shapely.MultiPolygon):
                    for poly in bub.polygon.geoms:
                        bubbles_refactored.append(
                            Bubble_v2(
                                polygon=poly, level=bub.level, is_hole=bub.is_hole
                            )
                        )
                else:
                    bubbles_refactored.append(bub)

            self.bubbles = bubbles_refactored
            self.__level2is_hole[level] = is_hole
            return True
        return False

    def intersection(
        *polygons: shapely.Polygon | shapely.MultiPolygon, grid_size
    ) -> shapely.Polygon | shapely.MultiPolygon:
        polygon = shapely.intersection_all([poly for poly in polygons])
        if polygon.area < grid_size:
            polygon = shapely.Polygon()
        else:
            if isinstance(polygon, shapely.Polygon | shapely.MultiPolygon):
                pass
            elif isinstance(polygon, shapely.GeometryCollection):
                polygon = shapely.MultiPolygon(
                    [p for p in polygon.geoms if p.area >= grid_size]
                )
            else:
                raise ValueError(f"polygon is of type {type(polygon)}")
        return polygon

    def difference(
        p: shapely.Polygon | shapely.MultiPolygon,
        q: shapely.Polygon | shapely.MultiPolygon,
        grid_size,
    ) -> shapely.Polygon | shapely.MultiPolygon:
        diff = p - q
        if diff.area < grid_size:
            diff = shapely.Polygon()
        else:
            if isinstance(diff, shapely.Polygon | shapely.MultiPolygon):
                pass
            elif isinstance(diff, shapely.GeometryCollection):
                diff = shapely.MultiPolygon(
                    [p for p in diff.geoms if p.area >= grid_size]
                )
            else:
                raise ValueError(f"diff is of type {type(diff)}")
        return diff

    def union(
        *polygons: shapely.Polygon | shapely.MultiPolygon, grid_size
    ) -> shapely.Polygon | shapely.MultiPolygon:
        polygon = shapely.union_all([poly for poly in polygons])
        if polygon.area < grid_size:
            polygon = shapely.Polygon()
        else:
            if isinstance(polygon, shapely.Polygon | shapely.MultiPolygon):
                pass
            else:
                raise ValueError(
                    f"polygon is of type {type(polygon)}, GeometryCollection should not appear here."
                )
        return polygon

    @staticmethod
    def __add_bubble_merge(bubbles, polygon, level, is_hole, grid_size):
        # Cut out higher level bubbles from new_bubble
        # TODO: bubble.level > new_bubble.level (>= ?)
        polygon = Topology_v2.difference(
            polygon,
            Topology_v2.union(
                *[
                    bubble.polygon
                    for bubble in bubbles
                    if bubble.level > new_bubble.level
                ],
                grid_size=grid_size,
            ),
            grid_size=grid_size,
        )
        # new_bubble is completely contained in higher level bubbles and is ignored
        if polygon.area < grid_size:
            return bubbles, False

        lower_bubbles = []
        disjoint_bubbles = []

        for bubble in bubbles:
            interserction = Topology_v2.intersection(
                polygon, bubble.polygon, grid_size=grid_size
            )
            if interserction.area < grid_size:
                disjoint_bubbles.append(bubble)
            else:
                if level > bubble.level:
                    diff = Topology_v2.difference(
                        bubble.polygon, new_bubble.polygon, grid_size=grid_size
                    )
                    if diff.area < grid_size:
                        # bubble is completely contained in new_bubble and is ignored
                        pass
                    else:
                        lower_bubbles.append(
                            Bubble_v2(
                                polygon=diff,
                                level=bubble.level,
                                is_hole=bubble.is_hole,
                            )
                        )
                        # TODO: Add intersection points to new_bubble
                        new_bubble.polygon = add_interserction_points(
                            new_bubble.polygon, bubble.polygon, grid_size=grid_size
                        )

                elif level == bubble.level:
                    # TODO: Check is_hole is equal
                    # TODO: Fix Multipolygon thingy (everywhere) or generalize the Bubbles class
                    new_bubble = Bubble_v2(
                        polygon=new_bubble.polygon.union(
                            bubble.polygon, grid_size=grid_size
                        ),
                        level=new_bubble.level,
                        is_hole=new_bubble.is_hole,
                    )
                else:
                    # Higher level bubbles are handeled in the beginning.
                    pass

        return disjoint_bubbles + lower_bubbles + [new_bubble], True

    def plot(self) -> Any:
        plot2(self.bubbles, self.level2is_hole)


class Topology_v3:
    DEFAULT_GRID_SIZE = 1e-15

    def __init__(
        self,
        mask: shapely.Polygon | shapely.MultiPolygon | None = None,
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        self.__grid_size = grid_size
        self.__topology = dict()
        self.__mask = mask

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def topology(self) -> dict[int, shapely.MultiPolygon]:
        return self.__topology

    @property
    def mask(self) -> shapely.Polygon | shapely.MultiPolygon | None:
        return self.__mask

    def add(
        self,
        multi_polygon: shapely.MultiPolygon,
        level: int,
        do_mask: bool = True,
        mask: shapely.MultiPolygon | shapely.Polygon = None,
        distances: dict[int, float] = dict(),
        transform: Callable[[shapely.MultiPolygon, Any], shapely.MultiPolygon] = None,
        transform_args: Any = None,
    ) -> bool:
        """Add new multipolygon to `self.topology`.

        Arguments:
            multi_polygon: shapely.MultiPolygon
                The multipolygon to be added to the topology.
            level: int
                Level of the multipolygon. Higher level polygons will be layered on top of lower level polygons.
                # TODO: Maybe rename to visibility_level.
            mask: shapely.MultiPolygon | shapely.Polygon
                If `None` self.mask is used, assuming it is not `None`. Only applies if `do_mask` is `True`.
            do_mask: bool
                If `True` the multipolygon is clipped with respect to the mask.
            distances: dict[int, float]
                Distances to other levels. If `None` all distances are set to 0.
            transform: Callable[[shapely.MultiPolygon, Any], shapely.MultiPolygon]
                Function that transforms the multipolygon before it is added to the topology.
                NOTE: Currently this only refers to periodicity.
            transform_args: Any
                Arguments for the transform function.

        Returns:
            `True` if the multipolygon was added to the topology, `False` otherwise.
        """
        # First transform the multipolygon if a transform function is given.
        if transform is not None:
            multi_polygon = (
                transform(multi_polygon, transform_args)
                if transform_args is not None
                else transform(multi_polygon)
            )
        # Apply mask if `do_mask` is `True` and a mask is given.
        if do_mask:
            if mask is None:
                if self.mask is None:
                    raise ValueError(
                        "No mask is given. Either set `do_mask` to `False` or set `self.mask`."
                    )
                else:
                    mask = self.mask
            # NOTE: multi_polygon might be a simple polygon after the intersection.
            multi_polygon = intersection(multi_polygon, mask, grid_size=self.grid_size)

        # Make sure the multipolygon is not too small
        if multi_polygon.area < self.grid_size:
            return False

        # Collision test
        if len(distances) > 0 and max(distances.values()) > self.grid_size:
            # Collision test is needed.
            for lvl, mp in self.topology.items():
                if lvl in distances:
                    if multi_polygon.distance(mp) < distances[lvl]:
                        return False

        # Add multipolygon to topology
        if level in self.topology:
            multi_polygon = union(
                self.topology[level], multi_polygon, grid_size=self.grid_size
            )
        if isinstance(multi_polygon, shapely.Polygon):
            # Make sure that the multipolygon is a MultiPolygon.
            print("WARNING: multi_polygon is a Polygon.")
            multi_polygon = shapely.MultiPolygon([multi_polygon])
        self.__topology[level] = multi_polygon
        return True

    def get_bubbles(self, cut_out_levels: list[int] = []) -> list[Bubble]:
        flatted_topology = self.topology
        levels = sorted(self.topology.keys())
        for i, lvl in enumerate(levels):
            # cut out higher level polygons from topology[lvl]
            for j in range(i + 1, len(levels)):
                diff = difference(
                    flatted_topology[lvl],
                    flatted_topology[levels[j]],
                    grid_size=self.grid_size,
                )
                if diff.area < self.grid_size:
                    flatted_topology[lvl] = shapely.Polygon()
                    break

                # Add intersection points to higher level polygons
                # TODO: Check if this works.
                flatted_topology[levels[j]] = intersection(
                    flatted_topology[levels[j]],
                    union(
                        flatted_topology[lvl],
                        flatted_topology[levels[j]],
                        grid_size=self.grid_size,
                    ),
                    grid_size=self.grid_size,
                )
                flatted_topology[lvl] = diff
                flatted_topology[lvl] = intersection(
                    flatted_topology[lvl],
                    union(
                        flatted_topology[levels[j]],
                        flatted_topology[lvl],
                        grid_size=self.grid_size,
                    ),
                    grid_size=self.grid_size,
                )

        bubbles = []
        for level, multipolygon in flatted_topology.items():
            if isinstance(multipolygon, shapely.MultiPolygon):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(
                                polygon=poly,
                                level=level,
                                is_hole=level in cut_out_levels,
                            )
                        )
            elif isinstance(multipolygon, shapely.Polygon):
                if multipolygon.area > self.grid_size:
                    bubbles.append(
                        Bubble(
                            polygon=multipolygon,
                            level=level,
                            is_hole=level in cut_out_levels,
                        )
                    )
            elif isinstance(multipolygon, shapely.GeometryCollection):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(
                                polygon=poly,
                                level=level,
                                is_hole=level in cut_out_levels,
                            )
                        )
            else:
                raise ValueError(f"multipolygon is of type {type(multipolygon)}")
        return bubbles

    def get_bubbles_old(self, cut_out_levels: list[int] = []) -> list[Bubble]:
        """Finalize the topology.

        Given polygons P, Q with P.level < Q.level, Q is cut out from P.
        Add intersection points to Q by doing:
        Q = intersection(Q, union(P, Q))
        => Get all intersection points by doing total_union = union(all polygons)
        iterate over all polygons and do:
            # either on polygon or on multipolygon
            polygon = intersection(polygon, total_union)
        total_union can be running_union.

        Arguments:
            cut_out_levels: list[int]
                Levels that should be cut out.
        """
        # bubbles = []
        # running_union = union(*self.topology.values(), grid_size=self.grid_size)

        # for level in sorted(self.topology.keys(), reverse=True):
        #     is_hole = level in cut_out_levels
        #     multi_polygon = intersection(self.topology[level], running_union, grid_size=self.grid_size)
        #     if isinstance(multi_polygon, shapely.Polygon):
        #         multi_polygon = shapely.MultiPolygon([multi_polygon])
        #     assert isinstance(multi_polygon, shapely.MultiPolygon)
        #     for polygon in multi_polygon.geoms:
        #         if polygon.area > self.grid_size:
        #             bubbles.append(Bubble(polygon=polygon, level=level, is_hole=is_hole))
        #     running_union = running_union - multi_polygon

        # return bubbles
        levels, multipolygons = zip(*sorted(self.topology.items()))
        multipolygons = list(multipolygons)

        n_multipolygons = len(multipolygons)
        lvl2intersection = dict()  # level to list of intersection points
        flatted_topology = dict()  # level to multipolygon

        for i in range(n_multipolygons):
            level = levels[i]
            remaining_multipolygon = multipolygons[i]
            is_empty = False
            # Cut out higher level bubbles from multipolygon.
            for j in range(i + 1, n_multipolygons):
                _diff = difference(
                    remaining_multipolygon, multipolygons[j], grid_size=self.grid_size
                )

                if _diff.area < self.grid_size:
                    is_empty = True
                    break
                # Add intersectin points
                multipolygons[j] = intersection(
                    multipolygons[j],
                    union(
                        remaining_multipolygon,
                        multipolygons[j],
                        grid_size=self.grid_size,
                    ),
                    grid_size=self.grid_size,
                )
                remaining_multipolygon = _diff
            if not is_empty:
                flatted_topology[level] = remaining_multipolygon
            else:
                continue

        bubbles = []
        is_hole = level in cut_out_levels
        for level, multipolygon in flatted_topology.items():
            if isinstance(multipolygon, shapely.MultiPolygon):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(polygon=poly, level=level, is_hole=is_hole)
                        )
            elif isinstance(multipolygon, shapely.Polygon):
                if multipolygon.area > self.grid_size:
                    bubbles.append(
                        Bubble(polygon=multipolygon, level=level, is_hole=is_hole)
                    )
            elif isinstance(multipolygon, shapely.GeometryCollection):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(polygon=poly, level=level, is_hole=is_hole)
                        )
            else:
                raise ValueError(f"multipolygon is of type {type(multipolygon)}")
        return bubbles

    @staticmethod
    def add_bubble(
        polygon, level, level_to_distance, topology, grid_size, mask=None
    ) -> tuple[dict[int, shapely.MultiPolygon], bool]:
        """
        Arguments:
            polygon: shapely.Polygon | shapely.MultiPolygon
            level: int
            level_to_distance: dict[int, float]
            topology: dict[int, shapely.MultiPolygon]
            mask: clip along the mask

        Assumptions:
            - Each Multipolygon in topology has no empty area and is a valid MultiPolygon.
        """
        # Distances must be either bigger than grid_size or 0 for merging.
        if not all(
            [dist >= grid_size or dist == 0 for dist in level_to_distance.values()]
        ):
            raise ValueError(
                "Distances must be either bigger than grid_size or 0 for merging."
            )

        # Clip the bubble with respect to the mask if given
        # TODO: Here the periodicity can be added
        if mask:
            polygon = Topology_v2.intersection(polygon, mask, grid_size=grid_size)

        # Make sure the bubble is not too small
        if polygon.area < grid_size:
            return topology, False

        # Collision test
        for lvl, dist in level_to_distance.items():
            if polygon.distance(topology[lvl]) < dist:
                return topology, False

        if level in topology:
            topology[level] = Topology_v2.union(
                topology[level], polygon, grid_size=grid_size
            )
        else:
            if isinstance(polygon, shapely.MultiPolygon):
                topology[level] = shapely.MultiPolygon([p for p in polygon.geoms])
            else:
                topology[level] = shapely.MultiPolygon([polygon])

        if isinstance(topology[level], shapely.Polygon):
            topology[level] = shapely.MultiPolygon([topology[level]])

        if not isinstance(topology[level], shapely.MultiPolygon):
            raise ValueError(f"topology[{level}] is of type {type(topology[level])}")

        return topology, True


def intersection(
    *polygons: shapely.Polygon | shapely.MultiPolygon, grid_size: float
) -> shapely.Polygon | shapely.MultiPolygon:
    """Return intersection of all polygons."""
    polygon = shapely.intersection_all([poly for poly in polygons])
    if polygon.area < grid_size:
        polygon = shapely.Polygon()
    else:
        if isinstance(polygon, shapely.Polygon | shapely.MultiPolygon):
            pass
        elif isinstance(polygon, shapely.GeometryCollection):
            polygon = shapely.MultiPolygon(
                [p for p in polygon.geoms if p.area >= grid_size]
            )
        else:
            raise ValueError(f"polygon is of type {type(polygon)}")
    return polygon


def union(
    *polygons: shapely.Polygon | shapely.MultiPolygon, grid_size: float
) -> shapely.Polygon | shapely.MultiPolygon:
    """Return the union of all polygons."""
    polygon = shapely.union_all([poly for poly in polygons])
    if polygon.area < grid_size:
        polygon = shapely.Polygon()
    else:
        if isinstance(polygon, shapely.Polygon | shapely.MultiPolygon):
            pass
        else:
            raise ValueError(
                f"polygon is of type {type(polygon)}, GeometryCollection should not appear here."
            )
    return polygon


def difference(
    p: shapely.Polygon | shapely.MultiPolygon,
    q: shapely.Polygon | shapely.MultiPolygon,
    grid_size: float,
) -> shapely.Polygon | shapely.MultiPolygon:
    """Return the difference p - q."""
    diff = p - q
    if diff.area < grid_size:
        diff = shapely.Polygon()
    else:
        if isinstance(diff, shapely.Polygon | shapely.MultiPolygon):
            pass
        elif isinstance(diff, shapely.GeometryCollection):
            diff = shapely.MultiPolygon([p for p in diff.geoms if p.area >= grid_size])
        else:
            raise ValueError(f"diff is of type {type(diff)}")
    return diff


def add_interserction_points(
    multi_poly: shapely.Polygon | shapely.MultiPolygon,
    poly: shapely.Polygon,
    grid_size: float,
) -> shapely.Polygon | shapely.MultiPolygon:
    """
    Add intersection to multi_poly.
    """
    return intersection(
        multi_poly, union(poly, multi_poly, grid_size=grid_size), grid_size=grid_size
    )


class Domain:
    DEFAULT_GRID_SIZE = 1e-15

    def __init__(
        self,
        reference_domain: shapely.Polygon | shapely.MultiPolygon,
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        self.__grid_size = grid_size
        self.__topology_constrained = dict()
        self.__ref_domain = reference_domain

        self.__topology_unconstrained = (
            dict()
        )  # level to unconstrained bubbles  used if stuff is removed

    def remove(self, level: int):
        pass


def is_valid(entity):
    assert isinstance(entity, shapely.Polygon) or isinstance(
        entity, shapely.MultiPolygon
    )


class Transform(object):
    def __call__(
        self, polygon: shapely.Polygon | shapely.MultiPolygon
    ) -> shapely.Polygon | shapely.MultiPolygon:
        pass


class PeriodicXY(Transform):
    def __init__(self, midpoint, width, height):
        pass


class PeriodicX(Transform):
    def __init__(self, midpoint, width, height):
        pass


class PeriodicY(Transform):
    def __init__(self, midpoint, width, height):
        pass


class Distance_Constraint(object):
    def set_level_distance(self, lvl_set_1, lvl_set_2, distance):
        pass

    def set_level_all(self, infos: dict()):
        pass

    def set_geometry_distance(self, level: shapely.Geometry, distance):
        pass

    def __call__(self, polygon):
        pass


"""distance_constraint = Distance_Constraint()
distance_constraint.set_level_distance(1,1,delta) 

  {1:  {1 : delta}}

distance_constraint.set_level_distance(1,[2,3,4], 2*delta) 
distance_constraint.set_geometry_distance(1, topology.domain.boundary, delta2)
topo.add(polygon, 1 , False, distance_constraint = distance_constraint)"""


class Topology_new:
    DEFAULT_GRID_SIZE = 1e-15

    def __init__(
        self,
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        self.__grid_size = grid_size
        self.__topology = dict()
        # self.__topology_unconstrained = [] #dict()

        self.__hole_levels = set()

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def topology(self) -> dict[int, shapely.Polygon | shapely.MultiPolygon]:
        return self.__topology

    @property
    def mask(self) -> shapely.Polygon | shapely.MultiPolygon | None:
        return self.__mask

    def __update_topology(self, polygon, level, is_hole, topology, merge=False):
        if is_hole:
            # if level not yet associated to holes...
            if level not in self.__hole_levels:
                # make sure it has not been treated as non level yet.
                if level in self.topology:
                    raise ValueError(
                        "Given level for hole was already associated with inclusion"
                    )
            self.__hole_levels.add(level)
        else:
            # make sure added level is not associated to holes
            if level in self.__hole_levels:
                raise ValueError(
                    "Given level was already associated with hole and cannot be used for an inclusion."
                )

        if not merge:
            if level not in topology:
                topology[level] = []
            self.topology[level].append(Bubble_v2(polygon, level, is_hole))
        else:
            if level not in topology:
                self.__topology[level] = polygon
            else:
                self.__topology[level] = shapely.union(
                    [self.__topology[level], polygon]
                )

    def __get_higher_levels(self, level):
        return [lvl for lvl in self.__topology.keys() if lvl > level]

    def add(
        self,
        polygon: shapely.Polygon | shapely.MultiPolygon,
        level: int,
        is_hole: bool = False,
        mask: shapely.Polygon | shapely.MultiPolygon | None = None,
        distance_constrained: Callable[
            [shapely.Polygon | shapely.MultiPolygon, Topology_new], bool
        ] = None,
        transform: Callable[
            [shapely.Polygon | shapely.MultiPolygon], shapely.MultiPolygon
        ] = None,
    ) -> bool:
        if transform is not None:
            polygon = transform(polygon)

        if mask is not None:
            polygon = intersection(polygon, mask, grid_size=self.grid_size)

        # self.__update_topology(polygon, level, is_hole, self.__topology_unconstrained)
        if distance_constrained is not None:
            if distance_constrained(polygon, self):
                return False

        # level hiearchy based constraints applied
        for higher_level in self.__get_higher_levels(level):
            assert len(self.__topology[level]) == 1
            subdomain = self.__topology[level]

            # polygon = difference(
            #         subdomain,
            #         polygon,
            #         grid_size=self.grid_size,
            #     )
            #     if diff.area < self.grid_size:
            #         flatted_topology[lvl] = shapely.Polygon()
            #         break

            #     # Add intersection points to higher level polygons
            #     # TODO: Check if this works.
            #     flatted_topology[levels[j]] = intersection(
            #         flatted_topology[levels[j]],
            #         union(
            #             flatted_topology[lvl],
            #             flatted_topology[levels[j]],
            #             grid_size=self.grid_size,
            #         ),
            #         grid_size=self.grid_size,
            #     )
            #     flatted_topology[lvl] = diff

        # Make sure the multipolygon is not too small
        if polygon.area < self.grid_size:
            return False

        # at this point the polygon must be polygon or multi polygon
        assert is_valid(polygon)
        self.__update_topology(polygon, level, is_hole, self.__topology, merge=True)

        return True

    def get_bubbles(self, cut_out_levels: list[int] = []) -> list[Bubble]:
        flatted_topology = self.topology
        levels = sorted(self.topology.keys())
        for i, lvl in enumerate(levels):
            # cut out higher level polygons from topology[lvl]
            for j in range(i + 1, len(levels)):
                diff = difference(
                    flatted_topology[lvl],
                    flatted_topology[levels[j]],
                    grid_size=self.grid_size,
                )
                if diff.area < self.grid_size:
                    flatted_topology[lvl] = shapely.Polygon()
                    break

                # Add intersection points to higher level polygons
                # TODO: Check if this works.
                flatted_topology[levels[j]] = intersection(
                    flatted_topology[levels[j]],
                    union(
                        flatted_topology[lvl],
                        flatted_topology[levels[j]],
                        grid_size=self.grid_size,
                    ),
                    grid_size=self.grid_size,
                )
                flatted_topology[lvl] = diff

        bubbles = []
        for level, multipolygon in flatted_topology.items():
            if isinstance(multipolygon, shapely.MultiPolygon):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(
                                polygon=poly,
                                level=level,
                                is_hole=level in cut_out_levels,
                            )
                        )
            elif isinstance(multipolygon, shapely.Polygon):
                if multipolygon.area > self.grid_size:
                    bubbles.append(
                        Bubble(
                            polygon=multipolygon,
                            level=level,
                            is_hole=level in cut_out_levels,
                        )
                    )
            elif isinstance(multipolygon, shapely.GeometryCollection):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(
                                polygon=poly,
                                level=level,
                                is_hole=level in cut_out_levels,
                            )
                        )
            else:
                raise ValueError(f"multipolygon is of type {type(multipolygon)}")
        return bubbles

    def get_bubbles_old(self, cut_out_levels: list[int] = []) -> list[Bubble]:
        """Finalize the topology.

        Given polygons P, Q with P.level < Q.level, Q is cut out from P.
        Add intersection points to Q by doing:
        Q = intersection(Q, union(P, Q))
        => Get all intersection points by doing total_union = union(all polygons)
        iterate over all polygons and do:
            # either on polygon or on multipolygon
            polygon = intersection(polygon, total_union)
        total_union can be running_union.

        Arguments:
            cut_out_levels: list[int]
                Levels that should be cut out.
        """
        # bubbles = []
        # running_union = union(*self.topology.values(), grid_size=self.grid_size)

        # for level in sorted(self.topology.keys(), reverse=True):
        #     is_hole = level in cut_out_levels
        #     multi_polygon = intersection(self.topology[level], running_union, grid_size=self.grid_size)
        #     if isinstance(multi_polygon, shapely.Polygon):
        #         multi_polygon = shapely.MultiPolygon([multi_polygon])
        #     assert isinstance(multi_polygon, shapely.MultiPolygon)
        #     for polygon in multi_polygon.geoms:
        #         if polygon.area > self.grid_size:
        #             bubbles.append(Bubble(polygon=polygon, level=level, is_hole=is_hole))
        #     running_union = running_union - multi_polygon

        # return bubbles
        levels, multipolygons = zip(*sorted(self.topology.items()))
        multipolygons = list(multipolygons)

        n_multipolygons = len(multipolygons)
        lvl2intersection = dict()  # level to list of intersection points
        flatted_topology = dict()  # level to multipolygon

        for i in range(n_multipolygons):
            level = levels[i]
            remaining_multipolygon = multipolygons[i]
            is_empty = False
            # Cut out higher level bubbles from multipolygon.
            for j in range(i + 1, n_multipolygons):
                _diff = difference(
                    remaining_multipolygon, multipolygons[j], grid_size=self.grid_size
                )

                if _diff.area < self.grid_size:
                    is_empty = True
                    break
                # Add intersectin points
                multipolygons[j] = intersection(
                    multipolygons[j],
                    union(
                        remaining_multipolygon,
                        multipolygons[j],
                        grid_size=self.grid_size,
                    ),
                    grid_size=self.grid_size,
                )
                remaining_multipolygon = _diff
            if not is_empty:
                flatted_topology[level] = remaining_multipolygon
            else:
                continue

        bubbles = []
        is_hole = level in cut_out_levels
        for level, multipolygon in flatted_topology.items():
            if isinstance(multipolygon, shapely.MultiPolygon):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(polygon=poly, level=level, is_hole=is_hole)
                        )
            elif isinstance(multipolygon, shapely.Polygon):
                if multipolygon.area > self.grid_size:
                    bubbles.append(
                        Bubble(polygon=multipolygon, level=level, is_hole=is_hole)
                    )
            elif isinstance(multipolygon, shapely.GeometryCollection):
                for poly in multipolygon.geoms:
                    if poly.area > self.grid_size:
                        bubbles.append(
                            Bubble(polygon=poly, level=level, is_hole=is_hole)
                        )
            else:
                raise ValueError(f"multipolygon is of type {type(multipolygon)}")
        return bubbles

    @staticmethod
    def add_bubble(
        polygon, level, level_to_distance, topology, grid_size, mask=None
    ) -> tuple[dict[int, shapely.MultiPolygon], bool]:
        """
        Arguments:
            polygon: shapely.Polygon | shapely.MultiPolygon
            level: int
            level_to_distance: dict[int, float]
            topology: dict[int, shapely.MultiPolygon]
            mask: clip along the mask

        Assumptions:
            - Each Multipolygon in topology has no empty area and is a valid MultiPolygon.
        """
        # Distances must be either bigger than grid_size or 0 for merging.
        if not all(
            [dist >= grid_size or dist == 0 for dist in level_to_distance.values()]
        ):
            raise ValueError(
                "Distances must be either bigger than grid_size or 0 for merging."
            )

        # Clip the bubble with respect to the mask if given
        # TODO: Here the periodicity can be added
        if mask:
            polygon = Topology_v2.intersection(polygon, mask, grid_size=grid_size)

        # Make sure the bubble is not too small
        if polygon.area < grid_size:
            return topology, False

        # Collision test
        for lvl, dist in level_to_distance.items():
            if polygon.distance(topology[lvl]) < dist:
                return topology, False

        if level in topology:
            topology[level] = Topology_v2.union(
                topology[level], polygon, grid_size=grid_size
            )
        else:
            if isinstance(polygon, shapely.MultiPolygon):
                topology[level] = shapely.MultiPolygon([p for p in polygon.geoms])
            else:
                topology[level] = shapely.MultiPolygon([polygon])

        if isinstance(topology[level], shapely.Polygon):
            topology[level] = shapely.MultiPolygon([topology[level]])

        if not isinstance(topology[level], shapely.MultiPolygon):
            raise ValueError(f"topology[{level}] is of type {type(topology[level])}")

        return topology, True


def union(
    *polygons: shapely.Geometry, grid_size
) -> shapely.Polygon | shapely.MultiPolygon:
    polygon = shapely.union_all([poly for poly in polygons])
    if polygon.area < grid_size:
        polygon = shapely.Polygon()
    else:
        if isinstance(polygon, shapely.Polygon | shapely.MultiPolygon):
            pass
        else:
            raise ValueError(
                f"polygon is of type {type(polygon)}, GeometryCollection should not appear here."
            )
    return polygon
