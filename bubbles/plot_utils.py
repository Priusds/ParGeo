"""Utility module, containing functions for plotting."""
from functools import cmp_to_key

import matplotlib.pyplot as plt
import shapely
from matplotlib.colors import Normalize


def plot(
    polygons: list[shapely.Polygon],
    polygon_levels: list[int],
    holes: set[int],
    domain: shapely.Polygon | shapely.MultiPolygon,
    diff_holes: bool = False,
) -> None:
    """
    Plot the topology.

    Args:
        polygons: list of polygons
        polygon_levels: list of levels for each polygon
        holes: set of levels that are holes
        domain: the domain
        diff_holes: If `True`, treat domain holes and holes with a level
            differently. If `False`, treat them the same, both filled white.
    """
    # Set of different levels.
    levels = set(polygon_levels)

    # Define all relevant colors.
    domain_boundary_color = "black"
    domain_color = "silver"
    hole_color = "white"
    interior_filling_color = "white"
    inter_hole_bound_color = "white"
    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            if lvl in holes and not diff_holes:
                lvl2cl[lvl] = hole_color
            else:
                lvl2cl[lvl] = plt.cm.cool(norm(lvl))

    # Make a legend.
    handles = []
    for lvl, color in lvl2cl.items():
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
    plt.legend(
        handles,
        [f"Level {lvl}" for lvl in lvl2cl.keys()],
        loc="upper left",
        bbox_to_anchor=(0.95, 1),
    )

    # Sort the polygons by inclusion.
    polygons, polygon_levels = zip(
        *sorted(
            zip(polygons, polygon_levels),
            key=cmp_to_key(polygon_compare),
            reverse=False,
        )
    )

    # Draw the reference domain's boundary.
    # TODO: Plot also the interior of holes
    for geo in domain.geoms:
        x, y = geo.exterior.xy
        plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
        for pp in geo.interiors:
            x, y = pp.xy
            plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
            plt.fill(x, y, color=interior_filling_color)

    for polygon, level in zip(polygons, polygon_levels):
        x, y = polygon.exterior.xy
        plt.fill(x, y, color=lvl2cl[level])
        if len(polygon.interiors) > 0:
            for pp in polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color=interior_filling_color)

        if level in holes:
            x, y = polygon.exterior.xy
            plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
            for sub_domain in domain.geoms:
                if polygon.intersects(sub_domain.exterior):
                    boundary_line = polygon.intersection(sub_domain.exterior)
                    if isinstance(boundary_line, shapely.MultiLineString):
                        for line in boundary_line.geoms:
                            x, y = line.xy
                            plt.plot(
                                x, y, "-", linewidth=2, color=inter_hole_bound_color
                            )
                    elif isinstance(boundary_line, shapely.LineString):
                        x, y = boundary_line.xy
                        plt.plot(x, y, "-", linewidth=2, color=inter_hole_bound_color)

            for interior in polygon.interiors:
                x, y = interior.xy
                plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)

    plt.show()


def polygon_compare(poly_1: shapely.Polygon, poly_2: shapely.Polygon) -> int:
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
            if shapely.Polygon(pp).contains(polygon_2):
                return -1
    if len(polygon_2.interiors) > 0:
        for pp in polygon_2.interiors:
            if shapely.Polygon(pp).contains(polygon_1):
                return 1

    return 0
