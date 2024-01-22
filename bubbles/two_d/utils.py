from functools import cmp_to_key
import shapely
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def plot(
    flattened_polygons: list[(shapely.Polygon, int)],
    holes: set[int],
    levels,
    domain,
    hole_color_mode="white",
):
    if len(holes) > 0:
        boundary_domain_poly = shapely.union_all(
            [
                poly
                for poly, lvl in flattened_polygons
                if (lvl not in holes and poly.intersects(domain.boundary))
            ]
        )
    sorted_bubbles = sorted(
        flattened_polygons, key=cmp_to_key(polygon_compare), reverse=False
    )

    # hole_boundary_color = "orange"
    domain_boundary_color = "black"
    domain_color = "silver"
    # ref_domain_hole_color = "white"

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            if lvl in holes and hole_color_mode == "white":
                lvl2cl[lvl] = "white"
            else:
                lvl2cl[lvl] = plt.cm.cool(norm(lvl))

    # Make a legend for the levels.
    handles = []
    for lvl, color in lvl2cl.items():
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color))
    plt.legend(
        handles,
        [f"Level {lvl}" for lvl in lvl2cl.keys()],
        loc="upper left",
        bbox_to_anchor=(0.95, 1),
    )

    # Draw the reference domain's boundary.
    # TODO: Plot also the interior of holes
    for geo in domain.geoms:
        x, y = geo.exterior.xy
        plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
        for pp in geo.interiors:
            x, y = pp.xy
            plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
            plt.fill(x, y, color="white")

    for polygon, level in sorted_bubbles:
        x, y = polygon.exterior.xy
        plt.fill(x, y, color=lvl2cl[level])
        if len(polygon.interiors) > 0:
            for pp in polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color="white")

        if level in holes:
            x, y = polygon.exterior.xy
            plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)
            for sub_domain in domain.geoms:
                if polygon.intersects(sub_domain.exterior):
                    boundary_line = polygon.intersection(sub_domain.exterior)
                    if isinstance(boundary_line, shapely.MultiLineString):
                        for line in boundary_line.geoms:
                            x, y = line.xy
                            plt.plot(x, y, "-", linewidth=2, color="white")
                    elif isinstance(boundary_line, shapely.LineString):
                        x, y = boundary_line.xy
                        plt.plot(x, y, "-", linewidth=2, color="white")
                    p = boundary_domain_poly.intersection(
                        boundary_line, grid_size=1e-15
                    )
                    p = (
                        p
                        if isinstance(p, shapely.MultiPoint)
                        else shapely.MultiPoint([p])
                    )
                    for pp in p.geoms:
                        x, y = pp.xy
                        plt.plot(
                            x,
                            y,
                            "s",
                            linewidth=2,
                            color=domain_boundary_color,
                            markersize=1,
                        )

            for interior in polygon.interiors:
                x, y = interior.xy
                plt.plot(x, y, "-", linewidth=2, color=domain_boundary_color)

    plt.show()


def polygon_compare(poly_1, poly_2):
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
