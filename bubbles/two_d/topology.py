"""Module for the topology of the two-dimensional domain.

Open Ideas:
- Create tree structure for the topology, useful for plotting and collision detection.
"""
from __future__ import annotations
from typing import Callable
import shapely
from bubbles.two_d.utils import plot


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
            holes: Levels that are holes. This can be changed at any time through `set_holes`.
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
            [shapely.Polygon | shapely.MultiPolygon, int, Topology], shapely.MultiPolygon
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
        if level <= 0:
            raise ValueError(f"Given level must be positive value, but given {level}.")

        # Apply transform if given
        if transform is not None:
            polygon = transform(polygon, level, self)

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

        # Make sure polygon intersects with the reference domain
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

                # TODO: Control intersection consistency.
                # assert( not isinstance(intersection, shapely.MultiPolygon) and not isinstance(intersection,shapely.Polygon))
                # TODO: intersection can be a polygon as well ! update logic here! in particular an empty polygon

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
                # TODO: avoid this: if statement
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

        # if level in self.__lvl2multipoly:
        #     if isinstance(self.__lvl2multipoly[level], shapely.GeometryCollection):
        #         raise ValueError(f"self.__lvl2multipoly[level] is a GeometryCollection")

        polygon = (
            polygon
            if not level in self.__lvl2multipoly
            else polygon.union(self.__lvl2multipoly[level], grid_size=self.grid_size)
        )

        # TODO: This case is really necessary, but should not be possible at this point !?
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

    def plot(self):
        """Plot the topology."""
        plot(
            self.flatten(),
            self.holes,
            self.levels,
            self.domain,
            hole_color_mode="white",
        )

    def flatten(self) -> list[(shapely.Polygon, int)]:
        flattened_polygons = []
        for level, running_polygon in self.__lvl2multipoly.items():
            if isinstance(running_polygon, shapely.MultiPolygon):
                for polygon in running_polygon.geoms:
                    flattened_polygons.append((polygon, level))
            else:
                flattened_polygons.append((running_polygon, level))
        return flattened_polygons
