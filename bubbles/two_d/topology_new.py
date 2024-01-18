from __future__ import annotations
from typing import Any, Callable
from pydantic import BaseModel
import shapely
from shapely.geometry import Polygon

from bubbles.two_d.hole import Rectangular

from shapely.affinity import translate as shapely_translate



class Topology:
    DEFAULT_GRID_SIZE = 1e-15

    def __init__(
        self,
        domain,
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        self.__grid_size = grid_size
        self.__topology = dict()
        self.__hole_levels = set()

        self.__reference_domain = domain

        # level = 0 associated with the domain
        self.__topology[0] = domain

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def topology(self) -> dict[int, shapely.Polygon | shapely.MultiPolygon]:
        return self.__topology
    
    @property
    def reference_domain(self) -> shapely.Polygon | shapely.MultiPolygon : 
        return self.__reference_domain

    # @property
    # def mask(self) -> shapely.Polygon | shapely.MultiPolygon | None:
    #     return self.__mask



    def __get_higher_levels(self, level):
        return [lvl for lvl in self.__topology.keys() if lvl > level]

    def add(
        self,
        polygon: shapely.Polygon | shapely.MultiPolygon,
        level: int,
        is_hole: bool = False,
        mask: shapely.Polygon | shapely.MultiPolygon | None = None,
        distance_constrained: Callable[
            [shapely.Polygon | shapely.MultiPolygon, Topology_new], bool
        ] = None,
        transform: Callable[
            [shapely.Polygon | shapely.MultiPolygon], shapely.MultiPolygon
        ] = None,
    ) -> bool:
        
        # check for valid level status
        self.__is_valid_level(level,is_hole)

        if transform is not None:
            polygon = transform(polygon)

        if mask is not None:
            polygon = shapely.intersection(polygon, self.reference_domain, grid_size=self.grid_size)
            # TODO: Intersection points needed ? 

        if distance_constrained is not None:
            if distance_constrained(polygon, self):
                return False

        # level hiearchy based constraints applied
        for higher_level in self.__get_higher_levels(level):
            assert len(self.__topology[level]) == 1
            
            # TODO: Intersect current polygon with self.__topology[higher_level]


        # Make sure the multipolygon is not too small
        if polygon.area < self.grid_size:
            return False

        # TODO: NOW is the correct time to check if intersection points have to be ADDED
        # TODO: Add intersection points possible to both polygon and self.__topology[higher_level]
        # e.g. intersection is a Linestring, single point, or polygon with additional points

        
        # at this point the polygon must be polygon or multi polygon
        assert is_valid(polygon)

        # The polygon can be merged/ added to the current topology. 
        self.__update_topology(polygon, level, is_hole, self.__topology, merge=True)

        return True
    



    def __is_valid_level(self, level, is_hole):

        if level == 0: 
            raise ValueError(
                        f"Given level = {level} is reserved for the underlying domain. Use level > 0 instead."
                    )

        if is_hole:
            # if level not yet associated to holes...
            if level not in self.__hole_levels:
                # make sure it has not been treated as non level yet.
                if level in self.topology:
                    raise ValueError(
                        "Given level for hole was already associated with inclusion"
                    )
            # mark level to be associated to hole_levels
            self.__hole_levels.add(level)
        else:
            # make sure added level is not associated to holes
            if level in self.__hole_levels:
                raise ValueError(
                    "Given level was already associated with hole and cannot be used for an inclusion."
                )


    def __update_topology(self, polygon, level, is_hole, topology, merge=False):
    
        if not merge:
            if level not in topology:
                topology[level] = []
            self.topology[level].append(Bubble_v2(polygon, level, is_hole))
        else:
            if level not in topology:
                self.__topology[level] = polygon
            else:
                self.__topology[level] = shapely.union(
                    [self.__topology[level], polygon]
                )