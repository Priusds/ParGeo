from __future__ import annotations
from typing import Any, Callable
from pydantic import BaseModel
import shapely
from shapely.geometry import Polygon

from bubbles.two_d.hole import Rectangular

from shapely.affinity import translate as shapely_translate


import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

class Subdomain(BaseModel):
    polygon: shapely.Polygon | shapely.MultiPolygon
    level: int
    is_hole: bool

    class Config:
        arbitrary_types_allowed = True

    # def __hash__(self) -> int:
    #     return self.polygon.__hash__()


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

        # clip with respect to the underlying reference domain
        polygon = shapely.intersection(polygon, self.reference_domain, grid_size=self.grid_size)
        # TODO: Intersection points needed ? needs to be done after the area test.

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
        # TODO: Add intersection points with on the level_0 reference domain if needed. 
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

    def plot(self): 
        pass


def plot(bubbles, level2is_hole, ref_domain) -> Any:
    """
    Plot the reference domain and the bubbles.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    bubbles = list(bubbles)
    bubbles.sort(key=lambda x: x.level)

    hole_boundary_color = "orange"
    ref_boundary_color = "black"
    ref_domain_color = "silver"
    ref_domain_hole_color = "white"

    # level2is_hole = self.level2is_hole TODO: Remove this

    # Filter out levels that are not represented by any bubble.
    # This is necessary as previous added bubbles might have been removed.
    present_levels = set([bubble.level for bubble in bubbles])
    level2is_hole = {
        lvl: is_hole for lvl, is_hole in level2is_hole.items() if lvl in present_levels
    }
    levels = set(level2is_hole.keys())

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = ref_domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            lvl2cl[lvl] = plt.cm.cool(norm(lvl))
    # Make a legend.
    handles = [
        Line2D([0], [0], label="Reference boundary", color=ref_boundary_color),
        Line2D([0], [0], label="Hole boundary", color=hole_boundary_color),
    ]
    for lvl, cl in lvl2cl.items():
        if lvl == 0:
            label_text = f"Domain lvl = {lvl}"
        else:
            label_text = (
                f"Hole lvl = {lvl}" if level2is_hole[lvl] else f"Inclusion lvl = {lvl}"
            )
        handles.append(mpatches.Patch(color=cl, label=label_text))
    plt.legend(handles=handles)

    # Plot the reference domain.
    x, y = ref_domain.polygon.exterior.xy
    plt.plot(x, y, color=ref_boundary_color)
    plt.fill(x, y, color=ref_domain_color)
    for pp in ref_domain.polygon.interiors:
        x, y = pp.xy
        plt.plot(x, y, color=ref_boundary_color)
        plt.fill(x, y, color=ref_domain_hole_color)

    # Plot the bubbles
    for bubble in bubbles:
        # Filter out reference domain
        if bubble.level != 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Re-Plot Bubbles without holes
    for bubble in bubbles:
        if len(bubble.polygon.interiors) == 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Show the plot
    plt.show()


def plot2(bubbles, level2is_hole) -> Any:
    """
    Plot the reference domain and the bubbles.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    bubbles.sort(key=lambda x: x.level)

    hole_boundary_color = "orange"
    ref_boundary_color = "black"
    ref_domain_color = "silver"
    ref_domain_hole_color = "white"

    # level2is_hole = self.level2is_hole TODO: Remove this

    # Filter out levels that are not represented by any bubble.
    # This is necessary as previous added bubbles might have been removed.
    present_levels = set([bubble.level for bubble in bubbles])
    level2is_hole = {
        lvl: is_hole for lvl, is_hole in level2is_hole.items() if lvl in present_levels
    }
    levels = set(level2is_hole.keys())

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = ref_domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            lvl2cl[lvl] = plt.cm.cool(norm(lvl))

    # Plot the bubbles
    for bubble in bubbles:
        x, y = bubble.polygon.exterior.xy
        plt.fill(x, y, color=lvl2cl[bubble.level])
        if len(bubble.polygon.interiors) > 0:
            for pp in bubble.polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color="white")

    # Re-Plot Bubbles without holes
    for bubble in bubbles:
        if len(bubble.polygon.interiors) == 0:
            x, y = bubble.polygon.exterior.xy
            plt.fill(x, y, color=lvl2cl[bubble.level])
            if bubble.is_hole:
                plt.plot(x, y, color=hole_boundary_color)

    # Show the plot
    plt.show()



def generate_topo():
    from bubbles.two_d.hole import Circle, Rectangular, Stellar, Ellipse

    from bubbles.two_d.topology import Topology as Topology_Working
    from bubbles.two_d.topology import Bubble
    import math
    """Generate a topology with a few edge cases."""
    # Make the domain
    d_1 = shapely.Polygon(
        Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    )
    d_2 = shapely.Polygon(
        Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)
    )

    dh_1 = shapely.Polygon(
        Circle(midpoint=(0.0, 0.75), radius=0.2).discretize_hole(refs=50)
    )

    dh_2 = shapely.Polygon(
        Circle(midpoint=(0.0, 0.0), radius=0.2).discretize_hole(refs=50)
    )
    dh_3 = shapely.Polygon(
        Ellipse(midpoint=(0.8, 0.0), axis=(0.05, 0.3), angle=0).discretize_hole(refs=50)
    )
    domain = d_1
    domain = domain.union(d_2)
    #domain = domain.difference(dh_1)
    #domain = domain.difference(dh_2)
    #domain = domain.difference(dh_3)
    topo = Topology_Working(domain)
    #topo.add(
    #    bubble=Bubble(polygon=domain, level=0, is_hole=False)
    #)
    
    # Add bubbles
    # c_1 = shapely.Polygon(
    #     Circle(midpoint=(-0.25, 0.0), radius=0.25).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_1, level=10, is_hole=True))
    
    
    # c_2 = shapely.Polygon(
    #     Circle(midpoint=(-0.1, 0.4), radius=0.3).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_2, level=11, is_hole=False))
    

    # c_3 = shapely.Polygon(
    #     Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_3, level=1, is_hole=False))
    

    # c_4 = shapely.Polygon(
    #     Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_4, level=1, is_hole=False))
    
    # # this guy is intersection with the domain
    # e_1 = shapely.Polygon(
    #     Ellipse(
    #         midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    #     ).discretize_hole(refs=200)
    # )
    # topo.add(Bubble(polygon=e_1, level=1, is_hole=False))
    

    # s_1 = shapely.Polygon(
    #     Stellar(midpoint=(-0.75, -0.6), radius=0.2).discretize_hole(refs=100)
    # )
    # topo.add(Bubble(polygon=s_1, level=1, is_hole=False))

    # s_2 = shapely.Polygon(
    #     Stellar(midpoint=(-0.5, -0.6), radius=0.25).discretize_hole(refs=100)
    # )
    # topo.add(Bubble(polygon=s_2, level=1, is_hole=False))
    

    # c_5 = shapely.Polygon(
    #     Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_5, level=3, is_hole=False))

    # c_6 = shapely.Polygon(
    #     Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_6, level=2, is_hole=False))
    
    c_7 = shapely.Polygon(
        Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_7, level=1, is_hole=False))
    
    is_hole = False
    level = 3
    e_2 = shapely.Polygon(
        Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
            refs=50
        )
    )
    topo.add(Bubble(polygon=e_2, level=level, is_hole=is_hole))

    e_3 = shapely.Polygon(
        Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
            refs=50
        )
    )
    topo.add(Bubble(polygon=e_3, level=level, is_hole=is_hole))

    e_4 = shapely.Polygon(
        Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
            refs=50
        )
    )
    topo.add(Bubble(polygon=e_4, level=level, is_hole=is_hole))

    # this guy encloses the loop !
    e_5 = shapely.Polygon(
        Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
            refs=50
        )
    )
    topo.add(Bubble(polygon=e_5, level=level, is_hole=is_hole))

    # c_8 = shapely.Polygon(
    #     Circle(midpoint=(-0.5, 0.4), radius=0.45).discretize_hole(refs=50)
    # )
    # topo.add(Bubble(polygon=c_8, level=3, is_hole=False))

    # c_9 = shapely.Polygon(
    #     Circle(midpoint=(0.0, 0.75), radius=0.01).discretize_hole(refs=50)
    # )

    # topo.add(Bubble(polygon=c_9, level=3, is_hole=False))

    return topo


def test_sorting():
        from functools import cmp_to_key
 
        def polygon_compare(poly_with_id_1, poly_with_id_2):
            polygon_1, id_1 = poly_with_id_1
            polygon_2, id_2 = poly_with_id_2

            if len(polygon_1.interiors) > 0 : 
                for pp in polygon_1.interiors :
                    if pp.contains(polygon_2): 
                        return -1
                            
            if len(polygon_2.interiors) > 0 : 
                for pp in polygon_2.interiors :
                    if pp.contains(polygon_1): 
                        return 1 

            return id_1 < id_2 ###

        compare_key = cmp_to_key(polygon_compare)

        def polygon_compare_simple(poly_1, poly_2):
            polygon_1, _ = poly_1
            polygon_2, _ = poly_2

            polygon_1 = shapely.Polygon(polygon_1.exterior)
            polygon_2 = shapely.Polygon(polygon_2.exterior)
            if polygon_1.contains(polygon_2):
                return 1
            elif polygon_2.contains(polygon_1):
                return -1
            # if len(polygon_1.interiors) > 0 : 
            #     for pp in polygon_1.interiors :
            #         if pp.contains(polygon_2): 
            #             return -1
            # if len(polygon_2.interiors) > 0 : 
            #     for pp in polygon_2.interiors :
            #         if pp.contains(polygon_1): 
            #             return 1 
            return 0


        return polygon_compare_simple
        


        



if __name__ == "__main__":
    from functools import cmp_to_key
    topo = generate_topo()

    # print(topo.bubbles)

    #topo.plot()
    #exit()

    flattened_polygons = [ (bubble.polygon, bubble.level) for bubble in topo.bubbles]

    polygon_compare_simple = test_sorting()





    
    sorted_bubbles = sorted(flattened_polygons, key = cmp_to_key(polygon_compare_simple), reverse=True)

    for k, A_t in enumerate(sorted_bubbles):
        for  l, B_t in enumerate(sorted_bubbles):
            if k > l : 
                A, lvlA = A_t
                B, lvlB = B_t
                print(f"{lvlA} : {lvlB} -> {polygon_compare_simple(A_t,B_t) }")
                print(f" A has interior : {len(A.interiors)}")
                print(f" B has interior : {len(B.interiors)}")
                print("================================================")


    #exit()

    #hole_boundary_color = "orange"
    #ref_boundary_color = "black"
    ref_domain_color = "silver"
    ref_domain_hole_color = "white"
    levels = topo.levels

    # Create a colormap for the levels.
    lvl2cl = dict()
    if 0 in levels:
        lvl2cl[0] = ref_domain_color
        levels.remove(0)
    if len(levels) > 0:
        norm = Normalize(vmin=min(levels), vmax=max(levels))
        for lvl in levels:
            lvl2cl[lvl] = plt.cm.cool(norm(lvl))

    for polygon, level in sorted_bubbles: 
        x, y = polygon.exterior.xy
        plt.fill(x, y, color=lvl2cl[level])
        if len(polygon.interiors) > 0:
            for pp in polygon.interiors:
                x, y = pp.xy
                plt.fill(x, y, color="white")
                
    plt.show()