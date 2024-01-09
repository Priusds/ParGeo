import shapely
import math

from bubbles.two_d.hole import Rectangular, Circle, Ellipse

from collections import namedtuple
from copy import deepcopy

Bubble = namedtuple("Bubble", ["polygon", "level", "is_hole"])

GRID_SIZE = 1e-15


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

    """

    def __init__(self, polygon: shapely.Polygon, grid_size=GRID_SIZE):
        # polygons is a list of non-intersecting polygons (0-measure intersection)
        # polygon = shapely.set_precision(polygon, grid_size)
        self.bubbles = {Bubble(polygon, level=0, is_hole=False)}

        self.ref_domain = Bubble(polygon, level=0, is_hole=False)
        self.grid_size = grid_size
        # TODO: Check if it is enough just to set the grid_size level for all new geometries,
        # and then leave the parameter out for the rest of the code.

        self.__level2is_hole = {0: False}

    @property
    def domain(self):
        domain_polygons = [p.polygon for p in self.bubbles if p.level == 0]
        return shapely.MultiPolygon(domain_polygons)

    @property
    def level2is_hole(self):
        return self.__level2is_hole

    @property
    def levels(self):
        return sorted(set([p.level for p in self.bubbles]))

    def add_bubble(self, Q: shapely.Polygon, level: int, is_hole: bool):
        if not level > 0:
            raise ValueError("Level must be positive integer.")
        if level in self.level2is_hole and self.level2is_hole[level] != is_hole:
            raise ValueError("Bubbles with same level must have same is_hole value.")
        # clip q with respect to the reference domain
        new_bubble = Q.intersection(self.ref_domain.polygon, grid_size=self.grid_size)
        if not isinstance(new_bubble, shapely.Polygon):
            raise ValueError(
                f"new_bubble is not a Polygon, found type {type(new_bubble)}."
            )

        # Find intersecting polygons and group by level relation
        lower_bubbles = set()
        higher_bubbles = set()
        equal_bubbles = set()
        for bubble in self.bubbles:
            if bubble.polygon.intersects(new_bubble):
                if level > bubble.level:
                    lower_bubbles.add(bubble)
                elif level < bubble.level:
                    higher_bubbles.add(bubble)
                elif level == bubble.level:
                    equal_bubbles.add(bubble)
                    if not bubble.is_hole == is_hole:
                        raise ValueError(
                            "Bubbles with same level must have same is_hole value."
                        )
                else:
                    raise ValueError(
                        "Given level is not an integer with order relation."
                    )
        if len(higher_bubbles) > 0:
            # Update the shape of new_bubble constrained by higher level bubbles
            # If new_bubble is completely contained within the union of higher level bubbles -> ignore q
            higher_polygons_union = shapely.union_all(
                [b.polygon for b in higher_bubbles]
            )
            if higher_polygons_union.contains(new_bubble):
                return
            new_bubble = new_bubble.difference(
                higher_polygons_union, grid_size=self.grid_size
            )
            # new_bubble might be a multipolygon
            if not (
                isinstance(new_bubble, shapely.Polygon)
                or isinstance(new_bubble, shapely.MultiPolygon)
            ):
                raise ValueError(f"new_bubble is of type {type(new_bubble)}")

        for lower_bubble in lower_bubbles:
            p_diff = lower_bubble.polygon.difference(
                new_bubble, grid_size=self.grid_size
            )
            if not (
                isinstance(p_diff, shapely.Polygon)
                or isinstance(p_diff, shapely.MultiPolygon)
            ):
                raise ValueError(f"p_diff is of type {type(p_diff)}")

            if p_diff.area > 0:
                if isinstance(p_diff, shapely.MultiPolygon):
                    for p_sub in p_diff.geoms:
                        assert isinstance(p_sub, shapely.Polygon)
                        self.bubbles.add(
                            Bubble(
                                p_sub,
                                level=lower_bubble.level,
                                is_hole=lower_bubble.is_hole,
                            )
                        )
                        # Add intersection points to new_bubble
                        new_bubble = new_bubble.intersection(
                            p_sub.union(new_bubble, grid_size=self.grid_size),
                            grid_size=self.grid_size,
                        )
                else:
                    # Add intersection points to new_bubble
                    new_bubble = new_bubble.intersection(
                        p_diff.union(new_bubble, grid_size=self.grid_size),
                        grid_size=self.grid_size,
                    )
                    self.bubbles.add(
                        Bubble(
                            p_diff,
                            level=lower_bubble.level,
                            is_hole=lower_bubble.is_hole,
                        )
                    )
            else:
                pass
                # # lower_bubble is covered by new_bubble
                # if lower_bubble.level == 0:
                #     # new_bubble completely contains the domain, this is not allowed
                #     return
                # else:
                #     # Remove lower_bubble from self.bubbles
                #     pass
            self.bubbles.remove(lower_bubble)

        if not (
            isinstance(new_bubble, shapely.Polygon)
            or isinstance(new_bubble, shapely.MultiPolygon)
        ):
            raise ValueError(f"new_bubble is of type {type(new_bubble)}")
        # Merge q_diff with equal level bubbles, the result still might be a multipolygon
        for lower_bubble in equal_bubbles:
            self.bubbles.remove(lower_bubble)

        q_final = shapely.union_all(
            [new_bubble] + [p.polygon for p in equal_bubbles], grid_size=self.grid_size
        )
        if not (
            isinstance(q_final, shapely.Polygon)
            or isinstance(q_final, shapely.MultiPolygon)
        ):
            raise ValueError(f"q_final is of type {type(q_final)}")

        if isinstance(q_final, shapely.MultiPolygon):
            for q_sub in q_final.geoms:
                self.bubbles.add(Bubble(q_sub, level=level, is_hole=is_hole))
        else:
            self.bubbles.add(Bubble(q_final, level=level, is_hole=is_hole))
        self.__level2is_hole[level] = is_hole

    def plot(self):
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

        # Plot the boundary of the reference domain.
        x, y = self.ref_domain.polygon.exterior.xy
        plt.plot(x, y, color=ref_boundary_color)
        for pp in self.ref_domain.polygon.interiors:
            x, y = pp.xy
            plt.plot(x, y, color=ref_boundary_color)

        # Plot the bubbles
        for bubble in bubbles:
            p, lvl, is_hole = bubble
            x, y = p.exterior.xy
            plt.fill(x, y, color=lvl2cl[lvl])
            if is_hole:
                plt.plot(x, y, color=hole_boundary_color)

        # Re-Plot Bubbles without holes
        for bubble in bubbles:
            p, lvl, is_hole = bubble
            if len(p.interiors) == 0:
                x, y = p.exterior.xy
                plt.fill(x, y, color=lvl2cl[lvl])
                if is_hole:
                    plt.plot(x, y, color=hole_boundary_color)

        # Show the plot
        plt.show()
