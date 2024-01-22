import shapely


class Constraints(object): 
    def __call__(self, polygon, level,  topology) -> True:
        pass
      

class ConvexityConstraint(Constraints): 
    pass


class DistanceConstraint(Constraints): 

    def __init__(self) -> None:
        self.distance_dict_poly = dict() # dict[(int,int) -> float]
        self.distance_dict_boundary = dict()
        self.geometry_distances = dict() # lvl -> list[(shapely.Geometry, float)]

        super().__init__()
    def _set_level_distance(self, lvl1, lvl2, dist, as_boundary_distance):

        if not as_boundary_distance: 
            self.distance_dict_poly[(lvl1, lvl2)] = dist
            if lvl1 != lvl2: 
                self.distance_dict_poly[(lvl2, lvl1)] = dist
        else: 
            self.distance_dict_boundary[(lvl1, lvl2)] = dist
            if lvl1 != lvl2: 
                self.distance_dict_boundary[(lvl2, lvl1)] = dist

    def set_level_distance(self, distance: float | list[float],
                                 level1 : int | list[int] | str = "any",
                                 level2 : int | list[int] | str = "any", 
                                 as_boundary_distance = False):
        
        levels_1 = None
        levels_2 = None
        distances= None
        if not isinstance(level1, str) and not isinstance(level2, str):
            levels_1 = [level1] if isinstance(level1,int) else level1 
            levels_2 = [level2] if isinstance(level2,int) else level2 
            distances= [distance] if isinstance(distance, float) else distance 

            if len(levels_1) != len(levels_2):  
                if len(levels_1) == 1: 
                    levels_1 = [levels_1[0]] * len(levels_2) 
                elif len(levels_2) == 1: 
                    levels_2 = [levels_2[0]] * len(levels_1)
                else: 
                    raise ValueError(f"Given level 1 = {level1} and level2 = {level2} not matching dimension, i.e. same length or one singleton.")

            if len(distances) != 1 and not len(distances) == len(levels_1):
                raise ValueError(f"Distances must have matching dimensions, i.e. a singleton or of matching length of level but got {distance}.")
            if len(distances) == 1: 
                distances = [distances[0]]*len(levels_1)
        else: 
            if level1 == "any" and not level2 == "any": 
                levels_2 = [level2] if isinstance(level2,int) else level2 
                levels_1 = ['any'] * len(levels_2)
                distances= [distance] if isinstance(distance, float) else distance 
            elif level2 == "any" and not level1 == "any":
                levels_1 = [level1] if isinstance(level1,int) else level1 
                levels_2 = ['any'] * len(levels_1)
                distances= [distance] if isinstance(distance, float) else distance 
            elif level1 == "any" and  level2 == "any":
                levels_1 = ['any']
                levels_2 = ['any']
                if not isinstance(distance, float) : 
                    raise ValueError("Distance parameter must be float in case of double presents of 'any' levels.")
                distances = [distance]

            if not (len(distances) == 1 or len(distances) == len(levels_1)): 
                raise ValueError("Distances has inconsistent dimension {len(distances)} for given input.")
            
            if len(distances) == 1 : 
                distances = [distances[0]]*len(levels_1)

        # levels_1, levels_2 and distances have same length now
        for lvl1, lvl2, d in zip(levels_1,levels_2, distances):
            self._set_level_distance(lvl1, lvl2, d, as_boundary_distance)

    

    def set_geometry_distance(self,  geometry: shapely.Geometry, distance : float | list[float], level: int | list[int] | str = "any"):
        levels    = []
        distances = []
        if not isinstance(level,str): 
            levels = [level] if isinstance(level,int) else level
            distances = [distance] if isinstance(distance,int) else distance
            if not len(levels) == 1  and len(distances) != len(levels): 
                raise ValueError(f"Distances must having matching dimension, i.e. 1 or same len as level given {len(levels)}.")
            if len(distances) == 1 : 
                distances = [distances[0]]*len(levels)
        else: 
            if not level == "any": 
                raise ValueError(f"Level as string must be 'any', but given as {level}.")
            distances = [distance] if isinstance(distance,float) else distance
            levels = ["any"]
            if not len(distances) == 1: 
                raise ValueError(f"Distances must be of length 1 for given level equals 'any' but given {len(distances)}")
            
        for lvl, d in zip(levels,distances):
            if level not in self.geometry_distances: 
                self.geometry_distances[level] = []
            self.geometry_distances[lvl].append((geometry, d))


    def __call__(self, polygon, level, topology):

        # TODO Can be computed once after each setter call
        semi_any_keys_poly = [keys[1] for keys in self.distance_dict_poly.keys() if 'any' == keys[0] and 'any' != keys[1]]
        semi_any_keys_boundary = [keys[1] for keys in self.distance_dict_boundary.keys() if 'any' == keys[0] and 'any' != keys[1]]

        # check if level distance constraint is violated to a polygon or boundary
        for level_multi_poly, multipoly in topology.lvl2multipoly.items():  

            if level_multi_poly == 0: 
                continue

            ## distance to given polygon distance constraints
            if ('any', 'any') in self.distance_dict_poly:
                if not distance_constraint(polygon, multipoly, self.distance_dict_poly[('any', 'any')], grid_size=topology.grid_size):
                    return False
                
            if ('any', level) in self.distance_dict_poly:
                if not distance_constraint(polygon, multipoly, self.distance_dict_poly[('any', level)], grid_size=topology.grid_size):
                    return False
            # level is part "any"  
            if level_multi_poly in semi_any_keys_poly: 
                if not distance_constraint(polygon, multipoly, self.distance_dict_poly[('any', level_multi_poly)], grid_size=topology.grid_size):
                    return False

            if (level_multi_poly, level) in self.distance_dict_poly:
                if not distance_constraint(polygon, multipoly, self.distance_dict_poly[(level_multi_poly, level)], grid_size=topology.grid_size):
                    return False
                
            ## distance to given boundary constraints
            if ('any', 'any') in self.distance_dict_boundary:
                if not distance_constraint(polygon, multipoly.boundary, self.distance_dict_boundary[('any', 'any')], grid_size=topology.grid_size):
                    return False
            if ('any', level) in self.distance_dict_boundary:
                if not distance_constraint(polygon, multipoly.boundary, self.distance_dict_boundary[('any', level)], grid_size=topology.grid_size):
                    return False
            # level is part "any"  
            if level_multi_poly in semi_any_keys_boundary: 
                if not distance_constraint(polygon, multipoly, self.distance_dict_boundary[('any', level_multi_poly)], grid_size=topology.grid_size):
                    return False
                
            if (level_multi_poly, level) in self.distance_dict_boundary:
                if not distance_constraint(polygon, multipoly.boundary, self.distance_dict_boundary[(level_multi_poly, level)], grid_size=topology.grid_size):
                    return False
        

        if level in self.geometry_distances: 
            for geo, dist in self.geometry_distances[level]:
                if not distance_constraint(polygon,geo, dist, grid_size=topology.grid_size):
                    return False
                
        if "any" in self.geometry_distances: 
            for geo, dist in self.geometry_distances["any"]:
                if not distance_constraint(polygon,geo, dist, grid_size=topology.grid_size):
                    return False

        
        # constraints are not violated
        return True

     


# TODO: Just add a collision: poly -> bool argument to the add function, it could for example
# check the number of refs.
def distance_constraint(
    poly: shapely.Polygon | shapely.MultiPolygon,
    geometry: shapely.Geometry,
    distance: float,
    grid_size: float = 1e-15,
):        
    """Check distance requirements between poly and geometry."""
    return poly.distance(geometry) > distance + grid_size