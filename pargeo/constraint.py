"""Constraint module for pargeo."""
from abc import ABC, abstractmethod
from itertools import product
from typing import Any, Dict

from shapely import Geometry, MultiPolygon, Polygon

from .topology import INITIAL_LEVEL, Topology


class Constraint(ABC):
    """Abstract constraint class."""

    @abstractmethod
    def __call__(
        self, polygon: Polygon | MultiPolygon, level: int, topology: Topology
    ) -> bool:
        """Check if constraint is violated by the polygon."""
        raise NotImplementedError


class ConvexityConstraint(Constraint):
    """Convexity constraint."""

    pass


class DistanceConstraint(Constraint):
    """Distance constraint.

    TODO: Contains quite some `# type: ignore` statements because shapely's
    `Geometry` class is not typed. Maybe there is a better way to do this.
    """

    def __init__(self) -> None:
        """Initialize the distance constraint.

        The distance constraint are stored in three dictionaries:
        - `__topo`: Distance constraints between topology multipolygons.
        - `__topo_bound`: Distance constraints between topology
            multipolygons boundaries.
        - `__geoms`: Distance constraints to shapely Geometry objects.
        """
        self.__topo: Dict[tuple[int | str, int | str], float] = dict()
        self.__topo_bound: Dict[tuple[int, int], float] = dict()
        self.__geoms: Dict[tuple[int, Geometry], float] = dict()

    def set_distance(
        self,
        distance: float | list[float],
        obj_1: int | str | list[int | str] = "any",
        obj_2: int | str | Geometry | list[int | str | Geometry] = "any",
        to_boundary=False,
    ):
        """Set variusous distance constraints.

        Distance constraints can be set between topology multipolygons with each
        other and between topology multipolygons and geometry objects. A geometry object
        is a shapely geometry object, e.g. shapely.Point, shapely.LineString, ect.

        Args:
            distance: Required distance between `obj_1` and `obj_2`.
                If a list of distances is given, the length of the list must
                either match the length of `obj_1` or `obj_2` or be single-valued.

            obj_1: The first object. Allowed values are:
                - int: The level of the topology multipolygon.
                - str: The string "any" to indicate any level.
                - list[int | str]: A list of the above values.

            obj_2: The second object. Allowed values are:
                - int: The level of the topology multipolygon.
                - str: The string "any" to indicate any level.
                - Geometry: A shapely Geometry object.
                - list[int | str | Geometry]: A list of the above values.

            to_boundary: If `True`, the distance is measured to the boundary.
                This is useful if inclusions are allowed.
        """
        # Make sure level_1 and level_2 are lists.
        if isinstance(obj_1, (int, str)):
            obj_1 = [obj_1]
        if isinstance(obj_2, (int, str, Geometry)):
            obj_2 = [obj_2]

        # Make sure distance is a list of distances.
        if isinstance(distance, float):
            distance = [distance]

        # Clarify types.
        obj_1: list[int | str] = obj_1  # type: ignore
        obj_2: list[int | str | Geometry] = obj_2  # type: ignore
        distance: list[float] = distance  # type: ignore

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
        for (lvl1, lvl2), dist in self.__topo.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Boundaries:")
        for (lvl1, lvl2), dist in self.__topo_bound.items():
            print(f"    {lvl1} - {lvl2}: {dist}")
        print("  Geometries:")
        for (lvl1, geom), dist in self.__geoms.items():
            print(f"    {lvl1} - {geom}: {dist}")

    def _set_distance(
        self,
        lvl1: int | str,
        lvl2: int | str | Geometry,
        dist: float,
        to_boundary: bool,
    ):
        """Write distance constraint to `__dist_dict_poly`, `__dist_dict_boundary` or
        `__dist_geoms`.

        For commodity reasons, the "distance matrix" is fully filled, a triangle
        matrix would be sufficient.
        """
        if isinstance(lvl2, Geometry):
            self.__geoms[(lvl1, lvl2)] = dist  # type: ignore
            return
        if not to_boundary:
            self.__topo[(lvl1, lvl2)] = dist  # type: ignore
            if lvl1 != lvl2:
                self.__topo[(lvl2, lvl1)] = dist  # type: ignore
        else:
            self.__topo_bound[(lvl1, lvl2)] = dist  # type: ignore
            if lvl1 != lvl2:
                self.__topo_bound[(lvl2, lvl1)] = dist  # type: ignore

    def __call__(
        self, polygon: Polygon | MultiPolygon, level: int, topology: Topology
    ) -> bool:
        """Compute the distance constraint."""
        # Get all encountered levels.
        topo_levels = topology.levels
        if level not in topo_levels:
            topo_levels.append(level)

        # MultiPolygon collision check, ignore `INITIAL_LEVEL`.
        dist_dict_poly = self.__reduce_distances(self.__topo, topo_levels, True)
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_poly.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            # Check collision level with lvl.
            if lvl in topology.lvl2multipoly and not check_distance_req(
                polygon,
                topology.lvl2multipoly[lvl],
                dist,
                grid_size=topology.grid_size,
            ):
                return False

        # MultiPolygon boundary collision check, do not ignore `INITIAL_LEVEL`.
        dist_dict_boundary = self.__reduce_distances(
            self.__topo_bound, topo_levels, False
        )
        lvls_dists = [
            (lvl2, dist)
            for (lvl1, lvl2), dist in dist_dict_boundary.items()
            if lvl1 == level
        ]
        for lvl, dist in lvls_dists:
            if lvl in topology.lvl2multipoly and not check_distance_req(
                polygon,
                topology.lvl2multipoly[lvl].boundary,
                dist,
                grid_size=topology.grid_size,
            ):
                return False

        # Geometry collision check, ignore `INITIAL_LEVEL`.
        dist_geoms = self.__reduce_distances(self.__geoms, topo_levels, True)
        geoms_dists = [
            (geom, dist) for (lvl1, geom), dist in dist_geoms.items() if lvl1 == level
        ]
        for geom, dist in geoms_dists:
            if not check_distance_req(
                polygon,
                geom,
                dist,
                grid_size=topology.grid_size,
            ):
                return False

        # All collision checks passed.
        return True

    @staticmethod
    def __reduce_distances(
        distance_dict: Dict[tuple[Any, Any], float],
        final_levels: list[int],
        ignore_initial_level: bool,
    ):
        """Resolve the `any` levels in the distance dictionary."""
        reduced_distance_dict: Dict[tuple[int, Any], float] = dict()
        if ignore_initial_level:

            def do_update(lvl_pair, d):
                return INITIAL_LEVEL not in lvl_pair and (
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


def check_distance_req(
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
