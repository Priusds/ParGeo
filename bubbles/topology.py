"""
Main module, use Topology to generate complex geometries.
"""
from __future__ import annotations

from typing import Callable, Optional

from shapely import GeometryCollection, MultiPolygon, Polygon

from bubbles.plot_utils import plot

DEFAULT_GRID_SIZE = 1e-15
INITIAL_LEVEL = 0


class Topology:
    """Class for the topology of the two-dimensional domain."""

    def __init__(
        self,
        domain: Polygon | MultiPolygon,
        holes: set[int] = set(),
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        """Initialize the topology.

        Args:
            domain: The domain.
            holes: Levels that are considered as holes.
                This can be changed at any time through `set_holes`.
            grid_size: The grid size.
        """
        if isinstance(domain, Polygon):
            domain = MultiPolygon([domain])
        self.__domain = domain
        self.__holes = holes  # Can be changed at any time.
        self.__grid_size = grid_size
        self.__lvl2multipoly = {INITIAL_LEVEL: domain}

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
    def domain(self) -> Polygon | MultiPolygon:
        return self.__domain

    @property
    def levels(self):
        return sorted(self.__lvl2multipoly.keys())

    @property
    def lvl2multipoly(self):
        return self.__lvl2multipoly

    def add(
        self,
        polygon: Polygon | MultiPolygon,
        level: int,
        transform: Optional[
            Callable[
                [Polygon | MultiPolygon, int, Topology],
                MultiPolygon,
            ]
        ] = None,
        constraint: Optional[
            Callable[[Polygon | MultiPolygon, int, Topology], bool]
        ] = None,
        extend_domain: bool = False,
    ) -> bool:
        """
        Add a polygon to the topology.

        Args:
            polygon: The polygon to add.
            level: The level of the polygon.
            transform: Is called to transform the polygon before adding it.
                If will call `transform(polygon, level, self)`.
            constraint: Is called to check if the polygon satisfies the constraints.
                The constraints are checked after the polygon is transformed.
                If the constraints are not satisfied, the polygon is not added.
                It will call `constraints(polygon, level, self)`.
            extend_domain: Flag to extend the domain if the polygon is outside of it.
        """
        if level < 0:
            raise ValueError("`level` must be non-negative.")

        # Apply transform if given
        if transform is not None:
            polygon = transform(polygon, level, self)

        # Clip the transformed polygon w.r.t. the domain
        if not extend_domain:
            polygon = polygon.intersection(self.domain, grid_size=self.grid_size)
        else:
            extended_domain = self.domain.union(polygon, grid_size=self.grid_size)
            if isinstance(extended_domain, Polygon):
                extended_domain = MultiPolygon([extended_domain])
            self.__domain = extended_domain

        # Make sure polygon is a MultiPolygon with non-zero area polygons.
        if isinstance(polygon, (GeometryCollection, MultiPolygon)):
            polygon = MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        # Make sure polygon intersects with the reference domain
        if polygon.area < self.grid_size:
            return False

        # TODO: Apply segmentize strategy

        # Check if the polygon satisfies the constraints
        if constraint is not None:
            if not constraint(polygon, level, self):
                return False

        # Update self.__lvl2multipoly, this is the main logic of the class.
        # Multipolygons with a higher level are cutted out from those with a lower
        # level, if the intersect. The level is a visibility ranking.
        new_lvl2multipoly = dict()
        for running_level, running_polygon in self.__lvl2multipoly.items():
            if not polygon.intersects(running_polygon):
                # running_polygon does not intersect with `polygon` so it can be
                # added directly to the new `__lvl2multipoly`.
                new_lvl2multipoly[running_level] = running_polygon
                continue

            # running_polygon intersects with `polygon`. Both might
            # need to be updated.
            if level != running_level:
                higher_polygon = polygon if level > running_level else running_polygon
                lower_polygon = polygon if level < running_level else running_polygon

                # Cut out the higher level polygon from the lower level polygon
                lower_polygon = lower_polygon.difference(
                    higher_polygon, grid_size=self.grid_size
                )

                # After shapely.difference, lower_polygon can be a GeometryCollection.
                if isinstance(lower_polygon, GeometryCollection):
                    lower_polygon = MultiPolygon(
                        [p for p in lower_polygon.geoms if p.area > self.grid_size]
                    )

                # If the lower polygon is approximately empty, it is not added.
                # If polygon = lower_polygon, then the polygon is not added.
                if lower_polygon.area < self.grid_size:
                    if level < running_level:
                        return False
                    continue

                # Add the intersection points, this is needed to create delauny
                # triangulations by GMSH.
                intersection_ = higher_polygon.intersection(
                    lower_polygon, grid_size=self.grid_size
                )

                # `intersection` can be a GeometryCollection, if so, iterate over the
                # geoms and add them one by one to the higher_polygon.
                if isinstance(intersection_, GeometryCollection):
                    for geom in intersection_.geoms:
                        # TODO: Check why geom can be a polygon, sometimes
                        # this happens, and the polygon has some approximate 0 area.

                        # Only add geoms with 0 area.
                        if geom.area == 0:
                            # Note geom.area == 0 is needed to remove polygons that,
                            # given our grid_size are not visible, *but* if a new
                            # polygon is added with some epsilon area, then this needs
                            # to be removed, which happens in the next if statement.
                            higher_polygon = higher_polygon.union(
                                geom, grid_size=self.grid_size
                            )
                            if isinstance(higher_polygon, GeometryCollection):
                                higher_polygon = MultiPolygon(
                                    [
                                        poly
                                        for poly in higher_polygon.geoms
                                        if poly.area > self.grid_size
                                    ]
                                )
                else:
                    higher_polygon = higher_polygon.union(
                        intersection_, grid_size=self.grid_size
                    )
                # TODO: Is it possible that the union in the upper if statement
                # is again a Polygon? If so, then we need to convert it.
                if isinstance(higher_polygon, Polygon):
                    higher_polygon = MultiPolygon([higher_polygon])

                # TODO: Is this necessary?
                if isinstance(higher_polygon, GeometryCollection):
                    higher_polygon = MultiPolygon(
                        [p for p in higher_polygon.geoms if p.area > self.grid_size]
                    )

                # Update the running level polygon and polygon.
                if running_level < level:
                    assert not isinstance(lower_polygon, GeometryCollection)
                    new_lvl2multipoly[running_level] = lower_polygon
                    polygon = higher_polygon  # polygon gets bigger

                else:  # running_level > level
                    assert not isinstance(higher_polygon, GeometryCollection)
                    new_lvl2multipoly[running_level] = higher_polygon
                    polygon = lower_polygon  # polygon gets smaller

        # Merge the polygon with the multipolygon of the same level.
        if level in self.__lvl2multipoly:
            polygon = polygon.union(
                self.__lvl2multipoly[level], grid_size=self.grid_size
            )
            # The shapely.union might create a GeometryCollection.
            # Appearently, polygon.union(multipolygon) can create
            # a GeometryCollection with a multiLineString.
            if isinstance(polygon, GeometryCollection):
                polygon = MultiPolygon(
                    [p for p in polygon.geoms if p.area > self.grid_size]
                )

        if isinstance(polygon, Polygon):
            polygon = MultiPolygon([polygon])

        if not isinstance(polygon, MultiPolygon):
            raise ValueError(f"`polygon` is not a MultiPolygon, but a {type(polygon)}.")

        new_lvl2multipoly[level] = polygon
        self.__lvl2multipoly = new_lvl2multipoly
        return True

    def plot(self) -> None:
        """Plot the topology."""
        plot(
            self.flatten(),
            self.holes,
            self.levels,
            self.domain,
            hole_color_mode="white",
        )

    def flatten(self) -> list[tuple[Polygon, int]]:
        """Return a list of polygons with their level."""
        flattened_polygons = []
        for level, running_polygon in self.__lvl2multipoly.items():
            if isinstance(running_polygon, MultiPolygon):
                for polygon in running_polygon.geoms:
                    flattened_polygons.append((polygon, level))
            else:
                flattened_polygons.append((running_polygon, level))
        return flattened_polygons
