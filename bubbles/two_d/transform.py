from typing import Any
import shapely
from shapely.affinity import translate as shapely_translate
from shapely.geometry import MultiPolygon, Polygon
import math


class Transform(object):
    """Base class for transforms."""

    def __call__(
        self, polygon, level, topology
    ) -> shapely.Polygon | shapely.MultiPolygon:
        pass

class Periodic(Transform):
    def __init__(self): #, levels: int | list[int] | str = "any"):
        #self.set_levels(levels)
        self.periodic_length_x = {}
        self.periodic_length_y = {}

    def set_periodicty(self, levels: int | list[int] | str = "any",
                             periodic_length_x : float | list[float] = float('inf'),
                             periodic_length_y : float | list[float] = float('inf')):
                    
        if isinstance(levels, str): 
            if not levels == "any": 
                raise ValueError(f"Levels specification as string must be 'any' but given {levels}")
            if not isinstance(periodic_length_x, float) or not isinstance(periodic_length_y, float): 
                raise ValueError(f"If given levels as 'any' periodic lengths must be given as a float")
            if not periodic_length_x > 0 or not periodic_length_y > 0: 
                raise ValueError(f"If given periodic length is supposed to be positive, but got {periodic_length_x} and {periodic_length_y}.")
                                 
            # any case overwrites all values ever set for fixed levels
            self.periodic_length_x = {"any": periodic_length_x}
            self.periodic_length_y = {"any": periodic_length_y}
            
        else: 
            levels_list = levels if isinstance(levels, list) else [levels]
            periodic_length_x_list = [periodic_length_x] * len(levels_list) if isinstance(periodic_length_x, float) else periodic_length_x
            periodic_length_y_list = [periodic_length_y] * len(levels_list) if isinstance(periodic_length_y, float) else periodic_length_y
            if not len(levels_list) == len(periodic_length_x_list):
                raise ValueError(f"Parameter periodic_length_x has no matching dimension, needed {len(levels_list)} but got {len(periodic_length_x_list)}")
            if not len(levels_list) == len(periodic_length_y_list):
                raise ValueError(f"Parameter periodic_length_y has no matching dimension, needed {len(levels_list)} but got {len(periodic_length_y_list)}")

            for lvl, Lx, Ly in zip(levels_list, periodic_length_x_list, periodic_length_y_list):
                if not Lx > 0 or not Ly > 0: 
                    raise ValueError(f"If given periodic length is supposed to be positive, but got {Lx} and {Ly}.")
                    
                self.periodic_length_x[lvl] = Lx 
                self.periodic_length_y[lvl] = Ly
                        

    def __call__(self, polygon, level, topology) -> shapely.Polygon | shapely.MultiPolygon:
        
        minx, miny, maxx, maxy = topology.domain.bounds

        poly_minx, poly_miny, poly_maxx, poly_maxy = polygon.bounds

        Lx = self.periodic_length_x["any"] if "any" in self.periodic_length_x else self.periodic_length_x[level]
        Ly = self.periodic_length_y["any"] if "any" in self.periodic_length_y else self.periodic_length_y[level]

        # poly_minx + kx_max * Lx <= maxx <=>   (maxx - poly_minx)/Lx >= kx_max 
        kx_max = math.floor((maxx - poly_minx)/Lx)  
        # poly_maxx + kx_min * Lx >= minx  <=>  (minx - poly_maxx)/Lx <= kx_min
        kx_min = math.ceil((minx - poly_maxx) / Lx)  
        # analog for y axis
        ky_max = math.floor((maxy - poly_miny)/Ly)  
        ky_min = math.ceil((miny - poly_maxy) / Ly)  

        # collect possible polygons for periodic rule
        periodic_polygons = [shapely.affinity.translate(polygon, xoff=kx * Lx, yoff= ky * Ly) for kx in range(kx_min, kx_max+1) for ky in range(ky_min, ky_max+1)]
        # and possible union them 
        periodic_polygons = shapely.union_all(periodic_polygons, 
                          grid_size=topology.grid_size)
        
        # TODO: If the periodic length matches with the boundary we possible need to add intersection points

        return periodic_polygons


# def clip_x(
#     polygon: shapely.Polygon | shapely.MultiPolygon,
#     x_min: float,
#     x_max: float,
#     grid_size: float = 1e-15,
# ):
#     """Clip the polygon along the x-axis."""
#     poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
#     if poly_x_min < x_min - grid_size:
#         poly_in = shapely.intersection(
#             polygon,
#             shapely.geometry.box(x_min, poly_y_min, poly_x_max, poly_y_max),
#             grid_size=grid_size,
#         )
#         poly_translated = shapely_translate(
#             shapely.intersection(
#                 polygon,
#                 shapely.geometry.box(poly_x_min, poly_y_min, x_min, poly_y_max),
#                 grid_size=grid_size,
#             ),
#             xoff=x_max - x_min,
#         )
#     elif poly_x_max > x_max + grid_size:
#         poly_in = shapely.intersection(
#             polygon,
#             shapely.geometry.box(poly_x_min, poly_y_min, x_max, poly_y_max),
#             grid_size=grid_size,
#         )
#         poly_translated = shapely_translate(
#             shapely.intersection(
#                 polygon,
#                 shapely.geometry.box(x_max, poly_y_min, poly_x_max, poly_y_max),
#                 grid_size=grid_size,
#             ),
#             xoff=x_min - x_max,
#         )
#     else:
#         assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
#         return polygon
#     if isinstance(poly_translated, shapely.GeometryCollection):
#         poly_translated = shapely.MultiPolygon(
#             [p for p in poly_translated.geoms if p.area > grid_size]
#         )
#     polygon = poly_in.union(poly_translated, grid_size=grid_size)
#     assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
#     return clip_x(polygon, x_min, x_max, grid_size=grid_size)


# def clip_y(
#     polygon: shapely.Polygon | shapely.MultiPolygon,
#     y_min: float,
#     y_max: float,
#     grid_size: float = 1e-15,
# ):
#     """Clip the polygon along the y-axis."""
#     poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
#     if poly_y_min < y_min - grid_size:
#         poly_in = shapely.intersection(
#             polygon,
#             shapely.geometry.box(poly_x_min, y_min, poly_x_max, poly_y_max),
#             grid_size=grid_size,
#         )
#         poly_translated = shapely_translate(
#             shapely.intersection(
#                 polygon,
#                 shapely.geometry.box(poly_x_min, poly_y_min, poly_y_max, y_min),
#                 grid_size=grid_size,
#             ),
#             yoff=y_max - y_min,
#         )
#     elif poly_y_max > y_max + grid_size:
#         poly_in = shapely.intersection(
#             polygon,
#             shapely.geometry.box(poly_x_min, poly_y_min, poly_x_max, y_max),
#             grid_size=grid_size,
#         )
#         poly_translated = shapely_translate(
#             shapely.intersection(
#                 polygon,
#                 shapely.geometry.box(poly_x_min, y_max, poly_x_max, poly_y_max),
#                 grid_size=grid_size,
#             ),
#             yoff=y_min - y_max,
#         )
#     else:
#         assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
#         return polygon
#     if isinstance(poly_translated, shapely.GeometryCollection):
#         poly_translated = shapely.MultiPolygon(
#             [p for p in poly_translated.geoms if p.area > grid_size]
#         )
#     polygon = poly_in.union(poly_translated, grid_size=grid_size)
#     assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
#     return clip_y(polygon, y_min, y_max, grid_size=grid_size)
