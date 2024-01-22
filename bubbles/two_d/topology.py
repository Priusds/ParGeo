from __future__ import annotations
from typing import Any, Callable
from pydantic import BaseModel
import shapely
from shapely.geometry import Polygon

from bubbles.two_d.hole import Rectangular

from shapely.affinity import translate as shapely_translate

from functools import cmp_to_key


import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


class Topology:
    DEFAULT_GRID_SIZE = 1e-15

    def __init__(
        self,
        domain,
        hole_levels: set[int] = set(),
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        self.__reference_domain = domain
        self.__grid_size = grid_size
        self.__topology = dict()
        self.__hole_levels = hole_levels

        # level = 0 associated with the domain
        self.__topology[0] = domain

    @property
    def hole_levels(self):
        return self.__hole_levels

    def set_hole_levels(self, hole_levels: set[int]):
        self.__hole_levels = hole_levels

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def topology(self) -> dict[int, shapely.Polygon | shapely.MultiPolygon]:
        return self.__topology

    @property
    def reference_domain(self) -> shapely.Polygon | shapely.MultiPolygon:
        return self.__reference_domain

    @property
    def levels(self):
        return sorted(self.topology.keys())

    def add(
        self,
        polygon: shapely.Polygon | shapely.MultiPolygon,
        level: int,
        distance_constrained: Callable[
            [shapely.Polygon | shapely.MultiPolygon, Topology], bool
        ] = None,
        transform: Callable[
            [shapely.Polygon | shapely.MultiPolygon], shapely.MultiPolygon
        ] = None,
    ) -> bool:
        # # check for valid level status
        # self.__is_valid_level(level,is_hole)
        # if is_hole:
        #     self.__hole_levels.add(level)

        if transform is not None:
            polygon = transform(polygon)

        # clip with respect to the underlying reference domain
        polygon = shapely.intersection(
            polygon, self.reference_domain, grid_size=self.grid_size
        )
        # Make sure polygon is not a GeometryCollection
        if isinstance(polygon, (shapely.GeometryCollection, shapely.MultiPolygon)):
            polygon = shapely.MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        # Make sure polygon intersects with the reference domain
        if polygon.area < self.grid_size:
            return False

        # TODO: Apply segmentize strategy

        if distance_constrained is not None:
            if distance_constrained(polygon, self):
                return False

        # apply geometrical constraints by level ordering
        updated_topo = dict()
        for running_level, running_polygon in self.topology.items():
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

                print("type intersection : ", type(intersection))
                print(intersection.area)
                print(intersection)
                #assert( not isinstance(intersection, shapely.MultiPolygon) and not isinstance(intersection,shapely.Polygon))
                #TODO: intersection can be a polygon as well ! update logic here! in particular an empty polygon
                

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
                #TODO: avoid this: if statement
                if isinstance(higher_polygon, shapely.GeometryCollection):
                    higher_polygon = shapely.MultiPolygon(
                        [p for p in higher_polygon.geoms if p.area > self.grid_size]
                    )

                # TODO: This might be unnecessary
                # lower_polygon = lower_polygon.union(
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


        # if level in self.topology:
        #     if isinstance(self.topology[level], shapely.GeometryCollection):
        #         raise ValueError(f"self.topology[level] is a GeometryCollection")

        polygon = (
            polygon
            if not level in self.topology
            else polygon.union(self.topology[level], grid_size=self.grid_size)
        )

        # TODO: This case is really necessary, but should not be possible at this point !?
        if isinstance(polygon, shapely.GeometryCollection):
            polygon = shapely.MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        updated_topo[level] = polygon

        # at this point the polygon must be polygon or multi polygon

        self.__topology = updated_topo
        return True

    def flatten(self):
        flattened_polygons = []
        for level, running_polygon in self.topology.items():
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
        holes_levels = self.hole_levels

        # hole_boundary_color = "orange"
        boundary_color = "black"
        ref_domain_color = "silver"
        # ref_domain_hole_color = "white"
        levels = self.levels

        reference_domain = self.reference_domain
        # Create a colormap for the levels.
        lvl2cl = dict()
        if 0 in levels:
            lvl2cl[0] = ref_domain_color
            levels.remove(0)
        if len(levels) > 0:
            norm = Normalize(vmin=min(levels), vmax=max(levels))
            for lvl in levels:
                lvl2cl[lvl] = plt.cm.cool(norm(lvl))

        # draw initial boundary
        x, y = reference_domain.exterior.xy
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
                if polygon.intersects(reference_domain.exterior):
                    boundary_line = polygon.intersection(reference_domain.exterior)
                    if isinstance(boundary_line, shapely.MultiLineString):
                        for k, line in enumerate(boundary_line.geoms):
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
