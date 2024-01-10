from typing import Any
from pydantic import BaseModel
import shapely


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

    def __add_bubble_disjoint(self, bubble: Bubble, bubble_distance: float, domain_distance: float) -> bool:
        """Add a bubble to the topology.
        
        The bubble is added only if it is disjoint to all other bubbles, having minimal distance `bubble_distance`.

        The bubble can intersect the domain boundary, if `domain_distance` is 0, otherwise it must be disjoint and have
        distance `domain_distance` to the domain boundary.
        """
        if not bubble.level > 0:
            raise ValueError("Level must be positive integer.")
        if bubble.level in self.level2is_hole and self.level2is_hole[bubble.level] != bubble.is_hole:
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
        self.bubbles.add(Bubble(
            polygon=polygon,
            level=bubble.level,
            is_hole=bubble.is_hole,
        ))
        self.__level2is_hole[bubble.level] = bubble.is_hole
        return True

    def add_bubble(self, bubble: Bubble, **kwargs) -> bool:
        """
        Add a bubble to the topology.
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
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D

        bubbles = list(self.bubbles)
        bubbles.sort(key=lambda x: x.level)

        hole_boundary_color = "orange"
        ref_boundary_color = "black"
        ref_domain_color = "silver"
        ref_domain_hole_color = "white"

        level2is_hole = self.level2is_hole
        # Filter out levels that are not represented by any bubble.
        # This is necessary as previous added bubbles might have been removed.
        present_levels = set([bubble.level for bubble in bubbles])
        level2is_hole = {
            lvl: is_hole
            for lvl, is_hole in level2is_hole.items()
            if lvl in present_levels
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
                    f"Hole lvl = {lvl}"
                    if level2is_hole[lvl]
                    else f"Inclusion lvl = {lvl}"
                )
            handles.append(mpatches.Patch(color=cl, label=label_text))
        plt.legend(handles=handles)

        # Plot the reference domain.
        x, y = self.ref_domain.polygon.exterior.xy
        plt.plot(x, y, color=ref_boundary_color)
        plt.fill(x, y, color=ref_domain_color)
        for pp in self.ref_domain.polygon.interiors:
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
