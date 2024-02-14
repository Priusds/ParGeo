"""
This module defines constraints for geometric objects in the pargeo library.

It includes an abstract base class `Constraint` that defines the interface for all constraints,
and `DistanceConstraint`.

`DistanceConstraint` allows setting distance constraints between different geometric objects.
These objects can be domain multipolygons, addressed by their level or shapely Geometry objects. 
The distances can be set between any two objects, or to the boundary of an object. 
The constraints are stored in three dictionaries, one for distances between domain multipolygons, one for distances to
the boundaries of domain multipolygons, and one for distances to shapely Geometry objects.


"""
from itertools import product
from typing import Any, Dict, Literal, overload

from pargeo.domain import Constraint, Domain
from pargeo.utils.typing_utils import (SHAPELY_GEOMETRIES, Level,
                                       ShapelyGeometry, SubDomain)

IntraConstraintDict = Dict[tuple[Level | Literal["any"], Level | Literal["any"]], float]

InterConstraintDict = Dict[tuple[Level | Literal["any"], ShapelyGeometry], float]


class DistanceConstraint(Constraint):
    """Distance constraint."""

    def __init__(self) -> None:
        """Initialize the distance constraint.

        The distance constraint are stored in three dictionaries:
        - `__polygon_constraints`: Distance constraints between the domain's subdomains.
        - `__boundary_constraints`: Distance constraints between domains' subdomain' boundaries.
        - `__geoms_constraints`: Distance constraints to arbitrary shapely Geometry objects.
        """
        self.__polygon_constraints: IntraConstraintDict = dict()

        self.__boundary_constraints: IntraConstraintDict = dict()

        self.__geoms_constraints: InterConstraintDict = dict()

    def set_distance(
        self,
        obj_1: Level | Literal["any"],
        obj_2: Level | Literal["any"] | ShapelyGeometry,
        distance: float,
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
        # Validate input.
        if not (isinstance(obj_1, int) or obj_1 == "any"):
            raise ValueError("`obj_1` must be a level or 'any'")

        if not (isinstance(obj_2, (int, SHAPELY_GEOMETRIES)) or obj_2 == "any"):
            raise ValueError("`obj_1` must be a level or 'any' or a shapely Geometry.")

        # Set distance constraints to shapely Geometry objects.
        if isinstance(obj_2, SHAPELY_GEOMETRIES):
            self.__geoms_constraints[(obj_1, obj_2)] = distance
            return

        # Set intra-domain distance constraints.
        # Distance constraint to the whole polygon.
        if not to_boundary:
            self.__polygon_constraints[(obj_1, obj_2)] = distance
            if obj_1 != obj_2:
                self.__polygon_constraints[(obj_2, obj_1)] = distance

        else:  # Distance constraint to the boundary of the polygon.
            self.__boundary_constraints[(obj_1, obj_2)] = distance
            if obj_1 != obj_2:
                self.__boundary_constraints[(obj_2, obj_1)] = distance

    def __call__(
        self, subdomain: SubDomain, level: Level, domain: Domain, **kwargs: Any
    ) -> bool:
        """Compute the distance constraint."""
        # Get all present levels of `domain`.
        levels = domain.levels
        grid_size = domain.grid_size
        background_level = domain.background_level
        if level not in levels:
            levels.append(level)

        # MultiPolygon collision check, ignores background.
        dist_dict_poly = self.__any_to_lvl(
            self.__polygon_constraints, levels, True, background_level
        )
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_poly.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            # Check collision level with lvl.
            if lvl in domain.level_to_subdomain:
                subdomain_2 = domain.level_to_subdomain[lvl]
                if subdomain.distance(subdomain_2) < dist - grid_size:
                    return False

        # MultiPolygon boundary collision check, DOESN'T ignore the backgroud.
        dist_dict_boundary = self.__any_to_lvl(
            self.__boundary_constraints, levels, False, background_level
        )
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_boundary.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            if lvl in domain.level_to_subdomain:
                subdomain_2 = domain.level_to_subdomain[lvl]
                if subdomain.distance(subdomain_2.boundary) < dist - grid_size:
                    return False

        # Geometry collision check, ignores the background.
        dist_geoms = self.__any_to_lvl(
            self.__geoms_constraints, levels, True, background_level
        )
        geoms_dists = [
            (geom, dist) for (lvl1, geom), dist in dist_geoms.items() if lvl1 == level
        ]
        for geom, dist in geoms_dists:
            if subdomain.distance(geom) < dist - grid_size:
                return False

        # All collision checks passed.
        return True

    def show(self):
        """Print the distance constraints."""
        print("Distance constraints:")
        print("Subdomains:")
        for (lvl1, lvl2), dist in self.__polygon_constraints.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Boundaries:")
        for (lvl1, lvl2), dist in self.__boundary_constraints.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Geometries:")
        for (lvl1, geom), dist in self.__geoms_constraints.items():
            print(f"    {lvl1} - {geom}: {dist}")

    @overload
    def __any_to_lvl(
        self,
        distance_dict: IntraConstraintDict,
        final_levels: list[Level],
        ignore_background: bool,
        background_level: Level,
    ) -> Dict[tuple[Level, Level], float]:
        ...

    @overload
    def __any_to_lvl(
        self,
        distance_dict: InterConstraintDict,
        final_levels: list[Level],
        ignore_background: bool,
        background_level: Level,
    ) -> Dict[tuple[Level, ShapelyGeometry], float]:
        ...

    def __any_to_lvl(
        self,
        distance_dict: IntraConstraintDict | InterConstraintDict,
        final_levels: list[Level],
        ignore_background: bool,
        background_level: Level,
    ) -> (
        Dict[tuple[Level, Level], float] | Dict[tuple[Level, ShapelyGeometry], float]
    ):  # TODO: Add type hints for return
        """Replace `any` with the actual levels."""
        reduced_distance_dict: Dict[tuple[Level, Level], float] | Dict[
            tuple[Level, ShapelyGeometry], float
        ] = dict()
        if ignore_background:

            def do_update(lvl1, lvl2, dist):
                return background_level not in (lvl1, lvl2) and (
                    lvl1,
                    lvl2 not in reduced_distance_dict
                    or reduced_distance_dict[lvl1, lvl2] < dist,
                )

        else:

            def do_update(lvl1, lvl2, dist):
                return (
                    lvl1,
                    lvl2 not in reduced_distance_dict
                    or reduced_distance_dict[lvl1, lvl2] < dist,
                )

        for (lvl_1, lvl_2), dist in distance_dict.items():
            match lvl_1, lvl_2:
                case "any", "any":
                    for lvl_1_, lvl_2_ in product(final_levels, repeat=2):
                        if do_update(lvl_1_, lvl_2_, dist):
                            reduced_distance_dict[lvl_1_, lvl_2_] = dist
                case "any", _:
                    for lvl in final_levels:
                        if do_update(lvl, lvl_2, dist):
                            reduced_distance_dict[(lvl, lvl_2)] = dist  # type: ignore
                case _, "any":
                    for lvl in final_levels:
                        if do_update(lvl_1, lvl, dist):
                            reduced_distance_dict[(lvl_1, lvl)] = dist  # type: ignore
                case _:
                    if do_update(lvl_1, lvl_2, dist):
                        reduced_distance_dict[(lvl_1, lvl_2)] = dist  # type: ignore

        return reduced_distance_dict
