"""Module for the topology of the two-dimensional domain.

Open Ideas:
- Create tree structure for the topology, useful for plotting and collision detection.
"""
from __future__ import annotations
from typing import Callable
import shapely
from shapely.affinity import translate as shapely_translate

from functools import cmp_to_key

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

DEFAULT_GRID_SIZE = 1e-15
DOMAIN_LEVEL = 0


class Topology:
    """Class for the topology of the two-dimensional domain."""

    def __init__(
        self,
        domain: shapely.Polygon | shapely.MultiPolygon,
        holes: set[int] = set(),
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        """Initialize the topology.

        Args:
            domain: The reference domain.
            holes: Levels that are holes. This can be changed at any time.
            grid_size: The grid size.
        """
        if isinstance(domain, shapely.Polygon):
            domain = shapely.MultiPolygon([domain])
        self.__domain = domain
        self.__holes = holes  # Can be changed at any time.
        self.__grid_size = grid_size
        self.__lvl2multipoly = {DOMAIN_LEVEL: domain}

    @property
    def holes(self):
        """Return the levels that are holes."""
        return self.__holes

    def set_holes(self, new_holes: set[int]):
        """Define the levels that are holes."""
        self.__holes = new_holes

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def domain(self) -> shapely.Polygon | shapely.MultiPolygon:
        return self.__domain

    @property
    def levels(self):
        return sorted(self.__lvl2multipoly.keys())

    def add(
        self,
        polygon: shapely.Polygon | shapely.MultiPolygon,
        level: int,
        domain_clip: bool = True,
        constraints: Callable[[shapely.Polygon | shapely.MultiPolygon], bool] = None,
        transform: Callable[
            [shapely.Polygon | shapely.MultiPolygon], shapely.MultiPolygon
        ] = None,
    ) -> bool:
        """
        Add a polygon to the topology.

        Args:
            polygon: The polygon to add.
            level: The level of the polygon.
            domain_clip: If True, clip the polygon w.r.t. the reference domain.
            constraints: A function that checks if the polygon satisfies some constraints.
                Example: The polygon must be inside the reference domain.
            transform: A function that transforms the polygon before adding it to the topology.
                Note: first the polygon is transformed and then the constraints are checked.
        """
        # Apply transform if given
        if transform is not None:
            polygon = transform(polygon)

        # Clip the transformed polygon w.r.t. the reference domain
        if domain_clip:
            polygon = shapely.intersection(
                polygon, self.domain, grid_size=self.grid_size
            )

        # Make sure polygon is not a GeometryCollection
        if isinstance(polygon, (shapely.GeometryCollection, shapely.MultiPolygon)):
            polygon = shapely.MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        # Make sure polygon is not empty
        if polygon.area < self.grid_size:
            return False

        # TODO: Apply segmentize strategy

        # Check if the polygon satisfies the constraints
        if constraints is not None:
            if not constraints(polygon):
                return False

        # apply geometrical constraints by level ordering
        updated_topo = dict()
        for running_level, running_polygon in self.__lvl2multipoly.items():
            if not polygon.intersects(running_polygon):
                updated_topo[running_level] = running_polygon
                continue

            if level != running_level:
                higher_polygon, higher_level = (
                    (polygon, level)
                    if level > running_level
                    else (running_polygon, running_level)
                )

                lower_polygon, lower_level = (
                    (polygon, level)
                    if level < running_level
                    else (running_polygon, running_level)
                )

                # subdomain can be a Point, LineString, MultiLineString, Polygon, Multipolygon or even GeometryCollection

                lower_polygon = lower_polygon.difference(
                    higher_polygon, grid_size=self.grid_size
                )

                if isinstance(lower_polygon, shapely.GeometryCollection):
                    lower_polygon = shapely.MultiPolygon(
                        [p for p in lower_polygon.geoms if p.area > self.grid_size]
                    )

                # in case of Point, Linestring, MultiLineString  or Polygon/MultiPolygon/GeometryCollection with area 0 dont add
                if lower_polygon.area < self.grid_size:
                    if level < running_level:
                        return False
                    continue

                # add the intersection points
                intersection = higher_polygon.intersection(
                    lower_polygon, grid_size=self.grid_size
                )

                if isinstance(intersection, shapely.GeometryCollection):
                    for geom in intersection.geoms:
                        if geom.area == 0:
                            # Note geom.area == 0 is needed to remove polygons that,
                            # given our grid_size are not visible, *but* if a new
                            # polygon is added with some epsilon area, then this needs to be
                            # removed, which happens in the next if statement.
                            higher_polygon = higher_polygon.union(
                                geom, grid_size=self.grid_size
                            )
                            if isinstance(higher_polygon, shapely.GeometryCollection):
                                higher_polygon = shapely.MultiPolygon(
                                    [
                                        p
                                        for p in higher_polygon.geoms
                                        if p.area > self.grid_size
                                    ]
                                )
                else:
                    higher_polygon = higher_polygon.union(
                        intersection, grid_size=self.grid_size
                    )
                if isinstance(higher_polygon, shapely.GeometryCollection):
                    higher_polygon = shapely.MultiPolygon(
                        [p for p in higher_polygon.geoms if p.area > self.grid_size]
                    )

                # TODO: This might be unnecessary
                # diff_polygon = diff_polygon.union(
                #     intersection, grid_size=self.grid_size
                # )

                # update the running level polygons
                if running_level == lower_level:
                    # update running level polygon shrinked by the higher level polygon that was added (with added intersection points)
                    updated_topo[running_level] = lower_polygon
                    # update the added polygon with added intersection points
                    polygon = higher_polygon

                elif running_level == higher_level:
                    # update running level polygon enriched by intersection points
                    updated_topo[running_level] = higher_polygon
                    # update the added polygon shrinked by the higher polygon (with added intersection points)
                    polygon = lower_polygon

        # Make sure the multipolygon is not too small
        # TODO: This might be unnecessary
        # if polygon.area < self.grid_size:
        #     return False

        if level in self.__lvl2multipoly:
            if isinstance(self.__lvl2multipoly[level], shapely.GeometryCollection):
                raise ValueError(f"self.__lvl2multipoly[level] is a GeometryCollection")

        polygon = (
            polygon
            if not level in self.__lvl2multipoly
            else polygon.union(self.__lvl2multipoly[level], grid_size=self.grid_size)
        )

        if isinstance(polygon, shapely.GeometryCollection):
            polygon = shapely.MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        if isinstance(polygon, shapely.Polygon):
            polygon = shapely.MultiPolygon([polygon])

        if not isinstance(polygon, shapely.MultiPolygon):
            raise ValueError(f"`polygon` is not a MultiPolygon, but a {type(polygon)}.")

        updated_topo[level] = polygon
        self.__lvl2multipoly = updated_topo
        return True

    def flatten(self):
        flattened_polygons = []
        for level, running_polygon in self.__lvl2multipoly.items():
            if isinstance(running_polygon, shapely.MultiPolygon):
                for polygon in running_polygon.geoms:
                    flattened_polygons.append((polygon, level))
            else:
                flattened_polygons.append((running_polygon, level))
        return flattened_polygons

    def plot(self):
        flattened_polygons = self.flatten()

        sorted_bubbles = sorted(
            flattened_polygons, key=cmp_to_key(polygon_compare), reverse=False
        )
        holes_levels = self.holes

        # hole_boundary_color = "orange"
        boundary_color = "black"
        ref_domain_color = "silver"
        # ref_domain_hole_color = "white"
        levels = self.levels

        domain = self.domain
        # Create a colormap for the levels.
        lvl2cl = dict()
        if 0 in levels:
            lvl2cl[0] = ref_domain_color
            levels.remove(0)
        if len(levels) > 0:
            norm = Normalize(vmin=min(levels), vmax=max(levels))
            for lvl in levels:
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
        # TODO: Fix AttributeError: 'MultiPolygon' object has no attribute 'exterior'
        x, y = domain.exterior.xy
        plt.plot(x, y, "-", linewidth=2, color=boundary_color)

        for polygon, level in sorted_bubbles:
            x, y = polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[level])
            if len(polygon.interiors) > 0:
                for pp in polygon.interiors:
                    x, y = pp.xy
                    plt.fill(x, y, color="white")

            if level in holes_levels:
                x, y = polygon.exterior.xy
                plt.plot(x, y, "-", linewidth=2, color=boundary_color)
                if polygon.intersects(domain.exterior):
                    boundary_line = polygon.intersection(domain.exterior)
                    if isinstance(boundary_line, shapely.MultiLineString):
                        for line in boundary_line.geoms:
                            x, y = line.xy
                            plt.plot(x, y, "-", linewidth=2, color="white")
                    elif isinstance(boundary_line, shapely.LineString):
                        x, y = boundary_line.xy
                        plt.plot(x, y, "-", linewidth=2, color="white")

        plt.show()


def polygon_compare(poly_1, poly_2):
    polygon_1, lvl1 = poly_1
    polygon_2, lvl2 = poly_2

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


def clip_x(
    polygon: shapely.Polygon | shapely.MultiPolygon,
    x_min: float,
    x_max: float,
    grid_size: float = 1e-15,
):
    """Clip the polygon along the x-axis."""
    poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
    if poly_x_min < x_min - grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(x_min, poly_y_min, poly_x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, poly_y_min, x_min, poly_y_max),
                grid_size=grid_size,
            ),
            xoff=x_max - x_min,
        )
    elif poly_x_max > x_max + grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, poly_y_min, x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(x_max, poly_y_min, poly_x_max, poly_y_max),
                grid_size=grid_size,
            ),
            xoff=x_min - x_max,
        )
    else:
        assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
        return polygon
    if isinstance(poly_translated, shapely.GeometryCollection):
        poly_translated = shapely.MultiPolygon(
            [p for p in poly_translated.geoms if p.area > grid_size]
        )
    polygon = poly_in.union(poly_translated, grid_size=grid_size)
    assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
    return clip_x(polygon, x_min, x_max, grid_size=grid_size)


def clip_y(
    polygon: shapely.Polygon | shapely.MultiPolygon,
    y_min: float,
    y_max: float,
    grid_size: float = 1e-15,
):
    """Clip the polygon along the y-axis."""
    poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
    if poly_y_min < y_min - grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, y_min, poly_x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, poly_y_min, poly_y_max, y_min),
                grid_size=grid_size,
            ),
            yoff=y_max - y_min,
        )
    elif poly_y_max > y_max + grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, poly_y_min, poly_x_max, y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, y_max, poly_x_max, poly_y_max),
                grid_size=grid_size,
            ),
            yoff=y_min - y_max,
        )
    else:
        assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
        return polygon
    if isinstance(poly_translated, shapely.GeometryCollection):
        poly_translated = shapely.MultiPolygon(
            [p for p in poly_translated.geoms if p.area > grid_size]
        )
    polygon = poly_in.union(poly_translated, grid_size=grid_size)
    assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
    return clip_y(polygon, y_min, y_max, grid_size=grid_size)


# TODO: Just add a collision: poly -> bool argument to the add function, it could for example
# check the number of refs.
def polygon_distance_req(
    poly1: shapely.Polygon | shapely.MultiPolygon,
    poly2: shapely.Polygon | shapely.MultiPolygon,
    distance: float,
    grid_size: float = 1e-15,
):
    """Check if two polygons are not too close."""
    return poly1.distance(poly2) > distance + grid_size


def boundary_distance_req(
    poly1: shapely.Polygon | shapely.MultiPolygon,
    poly2: shapely.Polygon | shapely.MultiPolygon,
    distance: float,
    grid_size: float = 1e-15,
):
    """Check if the boundary of two polygons are not too close."""
    return poly1.boundary.distance(poly2.boundary) > distance + grid_size
