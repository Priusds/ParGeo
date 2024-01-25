# from shapely.geometry import MultiPolygon, Polygon
import math
from typing import Any

import shapely
from shapely.affinity import translate as shapely_translate
from .topology import Topology

class Transform(object):
    """Base class for transforms."""

    def __call__(
        self, polygon, level, topology
    ) -> shapely.Polygon | shapely.MultiPolygon:
        pass


class Periodic(Transform):
    def __init__(self): 
        self.periodic_length_x = {}
        self.periodic_length_y = {}
        self.alpha = {}

    def set_periodicty(
        self,
        levels: int | list[int] | str = "any",
        periodic_length_x: float | list[float] = float("inf"),
        periodic_length_y: float | list[float] = float("inf"),
        alpha : float | list[float] = 0.
    ):

        if isinstance(levels, str):
            if not levels == "any":
                raise ValueError(
                    f"Levels specification as string must be 'any' but given {levels}"
                )
            if not isinstance(periodic_length_x, float) or not isinstance(
                periodic_length_y, float
            ):
                raise ValueError(
                    f"If given levels as 'any' periodic lengths must be given as a float"
                )
            if not periodic_length_x > 0 or not periodic_length_y > 0:
                raise ValueError(
                    f"If given periodic length is supposed to be positive, but got {periodic_length_x} and {periodic_length_y}."
                )

            # any case overwrites all values ever set for fixed levels
            self.periodic_length_x = {"any": periodic_length_x}
            self.periodic_length_y = {"any": periodic_length_y}
            self.alphas = {"any": alpha}

        else:
            levels_list = levels if isinstance(levels, list) else [levels]
            periodic_length_x_list = (
                [periodic_length_x] * len(levels_list)
                if isinstance(periodic_length_x, float)
                else periodic_length_x
            )
            periodic_length_y_list = (
                [periodic_length_y] * len(levels_list)
                if isinstance(periodic_length_y, float)
                else periodic_length_y
            )


            alpha_list = [alpha] * len(levels_list) if isinstance(alpha, float) else alpha

            if not len(levels_list) == len(periodic_length_x_list):
                raise ValueError(
                    f"Parameter periodic_length_x has no matching dimension, needed {len(levels_list)} but got {len(periodic_length_x_list)}"
                )
            if not len(levels_list) == len(periodic_length_y_list):
                raise ValueError(
                    f"Parameter periodic_length_y has no matching dimension, needed {len(levels_list)} but got {len(periodic_length_y_list)}"
                )

            if not len(levels_list) == len(alpha_list):
                raise ValueError(
                    f"Parameter alpha has no matching dimension, needed {len(levels_list)} but got {len(alpha_list)}"
                )

            for lvl, Lx, Ly, alph in zip(
                levels_list, periodic_length_x_list, periodic_length_y_list, alpha_list
            ):
                if not Lx > 0 or not Ly > 0:
                    raise ValueError(
                        f"If given periodic length is supposed to be positive, but got {Lx} and {Ly}."
                    )

                self.periodic_length_x[lvl] = Lx
                self.periodic_length_y[lvl] = Ly
                self.alphas[lvl] = alph

    def __call__(
        self, polygon, level, topology
    ) -> shapely.Polygon | shapely.MultiPolygon:

        alpha = (self.alphas["any"] 
        if "any" in self.alphas
            else self.alphas[level]
        )

        v = (math.cos(alpha), math.sin(alpha))
        w = (math.sin(alpha), -math.cos(alpha))

        Lv = (
            self.periodic_length_x["any"]
            if "any" in self.periodic_length_x
            else self.periodic_length_x[level]
        )
        Lw = (
            self.periodic_length_y["any"]
            if "any" in self.periodic_length_y
            else self.periodic_length_y[level]
        )
        
        # TODO: tightening this upper bound Kv and Kw if possible
        d = topology.domain.bounds[2] - topology.domain.bounds[0] + topology.domain.bounds[3] - topology.domain.bounds[1]
        Kv = math.ceil(d / Lv)
        Kw = math.ceil(d / Lw)
        
        v_dirs = [(v[0],v[1]), (-v[0],-v[1])]
        w_dirs = [(w[0],w[1]), (-w[0],-w[1])]

        periodic_polygons = shapely.union_all(
             [Repeat( v_dir, Lv, Kv, w_dir, Lw, Kw)(polygon, level, topology) for v_dir in v_dirs for w_dir in w_dirs], 
             grid_size=topology.grid_size)

        # TODO: If the periodic length matches with the boundary we possible need to add intersection points

        return periodic_polygons


class Repeat(Transform): 
        
    def __init__(self, 
            v_dir: tuple[float, float], 
            v_scalar: float,
            v_rep: int,
            w_dir: tuple[float, float]=None, 
            w_scalar: float = None, 
            w_rep: int = None,
            clip: bool = True,) -> Any:
            """Initialize a Repeat object.
        
            Args:
                v_dir: Direction vector. Will be normalized.
                v_scalar: Repeat distance.
                v_rep: Number of repeats.
                
                Note: Either all None or none of them.
                w_dir: Direction vector. Will be normalized.
                w_scalar: Repeat distance.
                w_rep: list[int]
            """
            if any((w_dir is None, w_scalar is None, w_rep is None)):
                if not all((w_dir is None, w_scalar is None, w_rep is None)):
                    raise ValueError("Only")
                
            self.v_dir = v_dir
            self.v_scalar = v_scalar
            self.v_rep = v_rep

            self.w_dir = w_dir
            self.w_scalar = w_scalar
            self.w_rep = w_rep

            self.v_data = (v_dir, v_scalar, v_rep)
            self.w_data = (w_dir, w_scalar, w_rep) if w_dir is not None else None
            
            self.clip = clip

        
    def __call__(self, polygon: shapely.Polygon, level: int, topo: Topology):
        """Repeated the polygon."""
        result_polygon = []

        bounds = polygon.bounds

        if self.w_data is None: 
            for i in range(self.v_rep + 1): 
                x_off = i * self.v_scalar * self.v_dir[0] 
                y_off = i * self.v_scalar * self.v_dir[1] 

                if self.clip : 
                    bounds_shifted = bounds[0] + x_off, bounds[1] + y_off, bounds[2] + x_off, bounds[3] + y_off
                    # fast sufficient check if the shifted polygon is outside domain
                    if not shapely.box(*bounds_shifted).intersects(topo.domain): 
                        continue

                poly_shifted = shapely.affinity.translate(polygon, xoff =x_off, yoff =  y_off)
                result_polygon.append(poly_shifted)
        else:
            for i in range(self.v_rep + 1):
                for j in range(self.w_rep + 1):
                    x_off = i * self.v_scalar * self.v_dir[0] + j * self.w_scalar * self.w_dir[0]
                    y_off = i * self.v_scalar * self.v_dir[1] + j * self.w_scalar * self.w_dir[1]

                    if self.clip : 
                        bounds_shifted = bounds[0] + x_off, bounds[1] + y_off, bounds[2] + x_off, bounds[3] + y_off
                        # fast sufficient check if the shifted polygon is outside domain
                        if not shapely.box(*bounds_shifted).intersects(topo.domain): 
                            continue
                    poly_shifted = shapely.affinity.translate(polygon, xoff =x_off, yoff =  y_off)                
                    result_polygon.append(poly_shifted)

        result = shapely.union_all(result_polygon)
        if not self.clip: 
            return result

        else: 
            intersection = result.intersection(topo.domain, grid_size = topo.grid_size)

            if isinstance(intersection, shapely.GeometryCollection):
                result = shapely.MultiPolygon([poly for poly in intersection if poly.area > topo.grid_size])

            return result


class Diffeomorphism(Transform):
    def __init__(self, mapping):
        """_summary_

        Args:
            mapping (_type_): _description_ Callable mapping IR^{2,batchsize} -> IR^2^{2,batchsize}
        """
        self.__mapping = mapping

    def __call__(
        self, polygon, level, topology
    ) -> shapely.Polygon | shapely.MultiPolygon:
        if isinstance(polygon, shapely.Polygon):
            return self.__map_polygon(polygon)
        elif isinstance(polygon, shapely.MultiPolygon):
            return shapely.MultiPolygon([self.__map_polygon(poly) for poly in polygon])

    def __map_polygon(self, polygon: shapely.Polygon):
        x, y = polygon.exterior.xy
        x_mapped, y_mapped = self.__mapping(x, y)

        shell = [(xx, yy) for xx, yy in zip(x_mapped, y_mapped)]

        holes = []
        if len(polygon.interiors) > 0:
            for hole in polygon.interiors:
                x, y = hole.exterior
                x_mapped, y_mapped = self.__mapping(x, y)
                hole_mapped = [(xx, yy) for xx, yy in zip(x_mapped, y_mapped)]
                holes.append(hole_mapped)
        return shapely.Polygon(shell, holes)
