"""
This module defines constraints for geometric objects in the pargeo library.

It includes an abstract base class `Constraint` that defines the interface for all constraints,
and two concrete constraint classes `ConvexityConstraint` and `DistanceConstraint`.

`DistanceConstraint` allows setting distance constraints between different geometric objects.
These objects can be domain multipolygons, addressed by their level or shapely Geometry objects. 
The distances can be set between any two objects, or to the boundary of an object. 
The constraints are stored in three dictionaries, one for distances between domain multipolygons, one for distances to
the boundaries of domain multipolygons, and one for distances to shapely Geometry objects.

`ConvexityConstraint` is currently without implementation.

"""
from itertools import product
from typing import Any, Dict, Literal

from shapely import Geometry, MultiPolygon, Polygon

from .domain import BACKGROUND_LEVEL, Constraint, Domain, Level, SubDomain


class DistanceConstraint(Constraint):
    """Distance constraint.

    TODO: Contains quite some `# type: ignore` statements because shapely's
    `Geometry` class is not typed. Maybe there is a better way to do this.
    """

    def __init__(self) -> None:
        """Initialize the distance constraint.

        The distance constraint are stored in three dictionaries:
        - `__polygon_constraints`: Distance constraints between the domain's subdomains.
        - `__boundary_constraints`: Distance constraints between domains' subdomain' boundaries.
        - `__geoms_constraints`: Distance constraints to arbitrary shapely Geometry objects.
        """
        self.__polygon_constraints: Dict[
            tuple[Level | Literal["any"], Level | Literal["any"]], float
        ] = dict()
        self.__boundary_constraints: Dict[tuple[Level, Level], float] = dict()
        self.__geoms_constraints: Dict[tuple[Level, Geometry], float] = dict()

    def set_distance(
        self,
        distance: float | list[float],
        obj_1: Level | Literal["any"] | list[Level | Literal["any"]],
        obj_2: Level
        | Literal["any"]
        | Geometry
        | list[Level | Literal["any"] | Geometry],
        to_boundary=False,
    ):
        """Set variusous distance constraints.

        Distance constraints can be set between domain multipolygons with each
        other and between domain multipolygons and geometry objects. A geometry object
        is a shapely geometry object, e.g. shapely.Point, shapely.LineString, ect.

        Args:
            distance: Required distance between `obj_1` and `obj_2`.
                If a list of distances is given, the length of the list must
                either match the length of `obj_1` or `obj_2` or be single-valued.

            obj_1: The first object. Allowed values are:
                - int: The level of the domain multipolygon.
                - str: The string "any" to indicate any level.
                - list[int | "any"]: A list of the above values.

            obj_2: The second object. Allowed values are:
                - int: The level of the domain multipolygon.
                - str: The string "any" to indicate any level.
                - Geometry: A shapely Geometry object.
                - list[int | str | Geometry]: A list of the above values.

            to_boundary: If `True`, the distance is measured to the boundary.
                This is useful if inclusions are allowed.
        """
        # Make sure level_1 and level_2 are lists.
        if isinstance(obj_1, (Level, str)):
            obj_1 = [obj_1]
        if isinstance(obj_2, (Level, str, Geometry)):
            obj_2 = [obj_2]

        # Make sure distance is a list of distances.
        if isinstance(distance, float):
            distance = [distance]

        # Verify input sizes.
        length_1 = len(obj_1)
        length_2 = len(obj_2)
        length_d = len(distance)
        if length_1 != length_2 and length_1 != 1 and length_2 != 1:
            raise ValueError
        if length_d not in {length_1, length_2} and length_d != 1:
            raise ValueError

        # Bring all to the same length to have a 1-to-1 relation between
        # level_1, level_2 and distance.
        max_length = max(length_1, length_2)
        if length_1 < max_length:
            obj_1 = obj_1 * max_length
        if length_2 < max_length:
            obj_2 = obj_2 * max_length
        if length_d < max_length:
            distance = distance * max_length

        # Set the distances.
        for lvl1, lvl2, d in zip(obj_1, obj_2, distance):
            self._set_distance(lvl1, lvl2, d, to_boundary)

    def show(self):
        """Print the distance constraints."""
        print("Distance constraints:")
        print("  Polygons:")
        for (lvl1, lvl2), dist in self.__polygon_constraints.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Boundaries:")
        for (lvl1, lvl2), dist in self.__boundary_constraints.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Geometries:")
        for (lvl1, geom), dist in self.__geoms_constraints.items():
            print(f"    {lvl1} - {geom}: {dist}")

    def _set_distance(
        self,
        lvl1: Level | Literal["any"],
        lvl2: Level | Literal["any"] | Geometry,
        dist: float,
        to_boundary: bool,
    ):
        """Write distance constraint to `__dist_dict_poly`, `__dist_dict_boundary` or
        `__dist_geoms`.

        For commodity reasons, the "distance matrix" is fully filled, a triangle
        matrix would be sufficient.
        """
        if isinstance(lvl2, Geometry):
            self.__geoms_constraints[(lvl1, lvl2)] = dist  # type: ignore
            return
        if not to_boundary:
            self.__polygon_constraints[(lvl1, lvl2)] = dist  # type: ignore
            if lvl1 != lvl2:
                self.__polygon_constraints[(lvl2, lvl1)] = dist  # type: ignore
        else:
            self.__boundary_constraints[(lvl1, lvl2)] = dist  # type: ignore
            if lvl1 != lvl2:
                self.__boundary_constraints[(lvl2, lvl1)] = dist  # type: ignore

    def __call__(
        self, subdomain: SubDomain, level: Level, domain: Domain, **kwargs: Any
    ) -> bool:
        """Compute the distance constraint."""
        # Get all encountered levels.
        levels = domain.levels
        if level not in levels:
            levels.append(level)

        # MultiPolygon collision check, ignore `BACKGROUND_LEVEL`.
        dist_dict_poly = self.__reduce_distances(
            self.__polygon_constraints, levels, True
        )
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_poly.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            # Check collision level with lvl.
            if (
                lvl in domain.level_to_subdomain
                and not DistanceConstraint.check_distance(
                    subdomain,
                    domain.level_to_subdomain[lvl],
                    dist,
                    grid_size=domain.grid_size,
                )
            ):
                return False

        # MultiPolygon boundary collision check, do not ignore `BACKGROUND_LEVEL`.
        dist_dict_boundary = self.__reduce_distances(
            self.__boundary_constraints, levels, False
        )
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_boundary.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            if (
                lvl in domain.level_to_subdomain
                and not DistanceConstraint.check_distance(
                    subdomain,
                    domain.level_to_subdomain[lvl].boundary,
                    dist,
                    grid_size=domain.grid_size,
                )
            ):
                return False

        # Geometry collision check, ignore `BACKGROUND_LEVEL`.
        dist_geoms = self.__reduce_distances(self.__geoms_constraints, levels, True)
        geoms_dists = [
            (geom, dist) for (lvl1, geom), dist in dist_geoms.items() if lvl1 == level
        ]
        for geom, dist in geoms_dists:
            if not DistanceConstraint.check_distance(
                subdomain,
                geom,
                dist,
                grid_size=domain.grid_size,
            ):
                return False

        # All collision checks passed.
        return True

    @staticmethod
    def __reduce_distances(
        distance_dict: Dict[tuple[Any, Any], float],
        final_levels: list[Level],
        ignore_background: bool,
    ):  # TODO: Add type hints for return
        """Resolve the `any` levels in the distance dictionary."""
        reduced_distance_dict: Dict[tuple[Level, Any], float] = dict()
        if ignore_background:

            def do_update(lvl_pair, d):
                return BACKGROUND_LEVEL not in lvl_pair and (
                    lvl_pair not in reduced_distance_dict
                    or reduced_distance_dict[lvl_pair] < d
                )

        else:

            def do_update(lvl_pair, d):
                return (
                    lvl_pair not in reduced_distance_dict
                    or reduced_distance_dict[lvl_pair] < d
                )

        for (lvl_1, lvl_2), dist in distance_dict.items():
            match lvl_1, lvl_2:
                case "any", "any":
                    for lvl_pair in product(final_levels, repeat=2):
                        if do_update(lvl_pair, dist):
                            reduced_distance_dict[lvl_pair] = dist  # type: ignore
                case "any", _:
                    for lvl in final_levels:
                        if do_update((lvl, lvl_2), dist):
                            reduced_distance_dict[(lvl, lvl_2)] = dist
                case _, "any":
                    for lvl in final_levels:
                        if do_update((lvl_1, lvl), dist):
                            reduced_distance_dict[(lvl_1, lvl)] = dist
                case _:
                    if do_update((lvl_1, lvl_2), dist):
                        reduced_distance_dict[(lvl_1, lvl_2)] = dist

        return reduced_distance_dict

    @staticmethod
    def check_distance(
        poly: Polygon | MultiPolygon,
        geometry: Geometry,
        distance: float,
        grid_size: float = 1e-15,
    ) -> bool:
        """Check distance requirements between poly and geometry.

        Returns:
            `False` if `poly` and `geometry` are too close.
        """
        return poly.distance(geometry) > distance + grid_size


class ConvexityConstraint(Constraint):
    """Convexity constraint."""

    pass
