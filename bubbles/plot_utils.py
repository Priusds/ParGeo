"""Utility module, containing functions for plotting."""
from functools import cmp_to_key

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from shapely import LineString, MultiLineString, MultiPolygon, Polygon


class DefaultColors:
    """Default colors for plotting."""

    boundary = "black"
    domain_boundary = "black"
    domain = "silver"
    hole = "white"
    interior_filling = "white"
    inter_hole_bound = "white"

    @staticmethod
    def get_color_map(levels, holes, color_holes=False):
        """Get a color map for the levels."""
        lvl2cl = dict()
        if 0 in levels:
            lvl2cl[0] = DefaultColors.domain
            levels.remove(0)
        if len(levels) > 0:
            norm = Normalize(vmin=min(levels), vmax=max(levels))
            for lvl in levels:
                if lvl in holes and not color_holes:
                    lvl2cl[lvl] = DefaultColors.hole
                else:
                    lvl2cl[lvl] = plt.cm.cool(norm(lvl))
        return lvl2cl


def make_legend(holes, colormap):
    """Make a legend."""
    handles = []
    descriptions = []
    for lvl, color in colormap.items():
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, edgecolor="black"))
        if lvl in holes:
            descriptions.append(f"{lvl}. level (hole)")
        else:
            descriptions.append(f"{lvl}. level")
    plt.legend(
        handles,
        descriptions,
        loc="upper left",
        bbox_to_anchor=(0.95, 1),
    ).set_draggable(True)


def plot_deprecated(
    polygons: list[Polygon],
    polygon_levels: list[int],
    holes: set[int],
    domain: Polygon | MultiPolygon,
    color_holes: bool = False,
) -> None:
    """
    Plot the topology.

    Args:
        polygons: list of polygons
        polygon_levels: list of levels for each polygon
        holes: set of levels that are holes
        domain: the domain
        color_holes: If `True` color added holes according to the level. Domain holes
            are colored white. If `False` color all holes white.
            TODO: Check if added polygons with holes are colored correctly.
    """
    # Set of different levels.
    levels = set(polygon_levels)

    color_map = DefaultColors.get_color_map(levels, holes, color_holes)

    make_legend(holes, color_map)

    # Sort the polygons by inclusion.
    polygons, polygon_levels = zip(  # type: ignore
        *sorted(
            zip(polygons, polygon_levels),
            key=cmp_to_key(__polygon_compare),
            reverse=False,
        )
    )

    # Draw the reference domain's boundary.
    # TODO: Plot also the interior of holes
    for geo in domain.geoms:
        x, y = geo.exterior.xy
        plt.plot(x, y, "-", linewidth=2, color=DefaultColors.domain_boundary)
        for pp in geo.interiors:
            x, y = pp.xy
            plt.plot(x, y, "-", linewidth=2, color=DefaultColors.domain_boundary)
            plt.fill(x, y, color=DefaultColors.interior_filling)

    for polygon, level in zip(polygons, polygon_levels):
        x, y = polygon.exterior.xy
        plt.fill(x, y, color=color_map[level])
        if len(polygon.interiors) > 0:
            for pp in polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color=DefaultColors.interior_filling)

        if level in holes:
            x, y = polygon.exterior.xy
            plt.plot(x, y, "-", linewidth=2, color=DefaultColors.domain_boundary)
            for sub_domain in domain.geoms:
                if polygon.intersects(sub_domain.exterior):
                    boundary_line = polygon.intersection(sub_domain.exterior)
                    if isinstance(boundary_line, MultiLineString):
                        for line in boundary_line.geoms:
                            x, y = line.xy
                            plt.plot(
                                x,
                                y,
                                "-",
                                linewidth=2,
                                color=DefaultColors.inter_hole_bound,
                            )
                    elif isinstance(boundary_line, LineString):
                        x, y = boundary_line.xy
                        plt.plot(
                            x, y, "-", linewidth=2, color=DefaultColors.inter_hole_bound
                        )

            for interior in polygon.interiors:
                x, y = interior.xy
                plt.plot(x, y, "-", linewidth=2, color=DefaultColors.domain_boundary)

    plt.show()


def __polygon_compare(poly_1: Polygon, poly_2: Polygon) -> int:
    """Compare two polygons."""
    polygon_1, _ = poly_1
    polygon_2, _ = poly_2

    # elements with interiors must always be before elements without interiors
    if len(polygon_1.interiors) > 0 and len(polygon_2.interiors) == 0:
        return -1
    if len(polygon_2.interiors) > 0 and len(polygon_1.interiors) == 0:
        return 1

    if len(polygon_1.interiors) > 0:
        for pp in polygon_1.interiors:
            if Polygon(pp).contains(polygon_2):
                return -1
    if len(polygon_2.interiors) > 0:
        for pp in polygon_2.interiors:
            if Polygon(pp).contains(polygon_1):
                return 1

    return 0
