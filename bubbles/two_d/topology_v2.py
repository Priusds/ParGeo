import shapely
import math

from bubbles.two_d.hole import Rectangular, Circle, Ellipse

from collections import namedtuple
from copy import deepcopy

Bubble = namedtuple("Bubble", ["polygon", "level", "is_hole"])


class Topology(object):
    """
    1. Add hole/inclusion q

    2. Find intersections I = p1,...,pn (elements from self.polygons)

    # update the shape of q constrained by higher level bubbles
    3. Find pj1, .., pjk from I such that all have higher level
    3.1 If q is completely contained within merge(pj1,...,pjk) -> ignore this hole
    3.2 modify q = q.difference(pj1, .., pjk) => q might be a multipolygon

    # update the lower level bubbles according to the (remaining) q
    4.  Find pj_1, .., pj_l from I such that all have lower level
        and modify pj_k = pj_k.difference(q) for all k = 1 .. l  => pj_k might be a multipolygon
    4.1 For each original pj_k that is fully contained in q: remove pj_k  ( maybe pj_k is empty set then ? )
    4.2 Else if pj_k is a Multipolygon then split it into multiple polygons and delete pj_k

    # At this point q might be a MultiPolygon
    5. Find pi1,...pim, subset of I, s.t. all have the same level as q
    5.1 Merge pi1,...pim and q. This is a (possible) multipolygon denoted as P.
    5.2 Delete pi1,...pim and add each subpolygon of P.
    5.3 If P is a multipolygon split into polygons.

    """

    def __init__(self, polygon: shapely.Polygon):
        # polygons is a list of non-intersecting polygons (0-measure intersection)
        self.bubbles = {Bubble(polygon, level=0, is_hole=False)}

        self.ref_domain = Bubble(polygon, level=0, is_hole=False)

        # HashMap:  level to list(polygons) map
        self.lvl2polygons = {}  # is final result for gmsh

    @property
    def domain(self):
        domain_polygons = [p.polygon for p in self.bubbles if p.level == 0]
        return shapely.MultiPolygon(domain_polygons)

    def set_domain(self, domain):
        pass

    def __add_discretized_void(self, discretized_hole, level):
        pass

    def __add_discretized_inclusion(self, discretized_hole, level, physical_group):
        pass

    def add_bubble(self, Q: shapely.Polygon, level: int, is_hole: bool):
        # clip q with respect to the reference domain
        q = Q.intersection(self.ref_domain.polygon)

        assert level > 0
        # Find intersecting polygons and group by level relation
        lower_bubbles = set()
        higher_bubbles = set()
        equal_bubbles = set()

        for p in self.bubbles:
            if p.polygon.intersects(q):
                if level > p.level:
                    lower_bubbles.add(p)
                elif level < p.level:
                    higher_bubbles.add(p)
                elif level == p.level:
                    equal_bubbles.add(p)
                else:
                    raise ValueError(
                        "Given level is not an integer with order relation."
                    )

        if len(higher_bubbles) > 0:
            print("I AM NEVER HERE")
            # Update the shape of q constrained by higher level bubbles
            # If q is completely contained within the union of higher level bubbles -> ignore q
            higher_polygons_union = shapely.union_all(
                [b.polygon for b in higher_bubbles]
            )
            if higher_polygons_union.contains(q):
                return
            # q_diff might be a multipolygon
            q = q.difference(higher_polygons_union)
            print(f"q_diff HL type: {type(q)}")
                        
        for p in lower_bubbles:

            p_diff = p.polygon.difference(q)
            if p_diff.area > 0:
                if isinstance(p_diff, shapely.MultiPolygon):
                    for p_sub in p_diff.geoms:
                        print(f"p_sub LL type: {type(p_sub)}")
                        self.bubbles.add(
                            Bubble(p_sub, level=p.level, is_hole=p.is_hole)
                        )
                else:
                    print(f"p_diff LL type: {type(p_diff)}")
                    self.bubbles.add(Bubble(p_diff, level=p.level, is_hole=p.is_hole))
            else:
                # p is covered by q_diff
                if p.level == 0:
                    # bubble completely contains the domain, this is not allowed
                    return
                else:
                    # Remove p from self.bubbles
                    pass
            self.bubbles.remove(p)

        # Clip q_diff with main domain
        #q_diff = q.intersection(self.domain)
        #print(f"q_diff intersection domain type: {type(q)}")
        # Merge q_diff with equal level bubbles, the result still might be a multipolygon
        q_final = shapely.union_all([q] + [p.polygon for p in equal_bubbles])
        #print(f"q_final domain type: {type(q_final)}")

        q_bubble = Bubble(q_final, level=level, is_hole=is_hole)

        # NOTE: currently the is_hole information is not used, this makes the implementation easier
        # but leads to redundant bubbles and more shapely operations.
        # IDEA: add an IF statement that only adds the new bubble if it is not a holes.
        self.bubbles.add(q_bubble)

    def add_void(self, discretize_hole, level):
        # assert hole has discretize method
        pass

    def add_inclusion(self, discretize_hole, level):
        # assert hole has discretize method
        pass

    def set_physical_group(self, level_to_physical_group_map):
        pass

    def plot_setup(self): 
        lvl2cl = {0: "blue", 1: "brown", 2: "gray", 3: "yellow", 10: "white", 11: "purple"}
        P = sorted(self.bubbles, key=lambda p: p.level)
        for bubble in P:
            p, lvl, is_hole = bubble 
            print(f"Plot bubble with level = {lvl}")
            x, y = p.exterior.xy  # gold
            if is_hole:
                plt.fill(x, y, color=lvl2cl[lvl])
                plt.plot(x, y, color='black')
            else:
                plt.fill(x, y, color=lvl2cl[lvl])

        
        x,y = self.ref_domain.polygon.exterior.xy
        plt.plot(x, y, color='black')
        for pp in self.ref_domain.polygon.interiors : 
            x, y = pp.xy  # gol
            plt.fill(x, y, color='white') 
            plt.plot(x, y, color='black')


        plt.show()


def test():
    lvl2cl = {0: "blue", 1: "brown", 2: "gray", 3: "yellow", 10: "white"}

    topo = Topology()

    R = Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    C = Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)
    topo.add_bubble(R, level=0, is_hole=False)
    topo.add_bubble(C, level=0, is_hole=False)

    C2 = Circle(midpoint=(-0.25, 0), radius=0.25).discretize_hole(refs=50)
    topo.add_bubble(C2, level=10, is_hole=True)

    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)

    topo.add_bubble(C3, level=1, is_hole=False)
    topo.add_bubble(C4, level=1, is_hole=False)

    # this guy is intersection with the domain
    E1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize_hole(refs=50)
    # this guy is just interior
    E2 = Ellipse(midpoint=(-0.55, -0.6), axis=(0.25, 0.1), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(E1, level=1, is_hole=False)
    topo.add_bubble(E2, level=1, is_hole=False)

    C5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    C6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    topo.add_bubble(C6, level=4, is_hole=False)
    topo.add_bubble(C5, level=3, is_hole=False)

    C = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    topo.add_bubble(C, level=1, is_hole=False)

    E = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(E, level=2, is_hole=False)

    E = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(E, level=2, is_hole=False)

    E = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(E, level=2, is_hole=False)

    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(E, level=2, is_hole=False)

    P = sorted(topo.polygons, key=lambda p: p.level)
    # check order

    for p in P:
        x, y = p.exterior.xy  # gold
        if p.is_hole:
            plt.fill(x, y, color=lvl2cl[p.level])
        else:
            plt.fill(x, y, color=lvl2cl[p.level])

    plt.show()


def test2():
    import matplotlib.pyplot as plt

    # Domain
    R = Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    C = Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)
    R_shapely = shapely.Polygon(R)
    C_shapely = shapely.Polygon(C)
    domain = R_shapely.union(C_shapely)

    C2 = Circle(midpoint=(-0.25, 0), radius=0.25).discretize_hole(refs=50)
    p2 = shapely.Polygon(C2)

    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    p3 = shapely.Polygon(C3)
    p4 = shapely.Polygon(C4)

    E1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize_hole(refs=50)
    p5 = shapely.Polygon(E1)

    E2 = Ellipse(midpoint=(-0.55, -0.6), axis=(0.25, 0.1), angle=0.0).discretize_hole(
        refs=50
    )
    p6 = shapely.Polygon(E2)

    C5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    C6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    p7 = shapely.Polygon(C5)
    p8 = shapely.Polygon(C6)

    domain = domain.difference(p2)

    x, y = domain.exterior.xy
    plt.fill(x, y, color="blue")

    x, y = domain.interiors[0].xy
    plt.fill(x, y, color="white")

    inclusion = p3.union(p4)
    x, y = inclusion.exterior.xy
    plt.fill(x, y, color="brown")

    p5 = p5.intersection(domain)
    x, y = p5.exterior.xy
    plt.fill(x, y, color="brown")

    x, y = p6.exterior.xy
    plt.fill(x, y, color="brown")

    x, y = p8.exterior.xy  # gold
    plt.fill(x, y, color="yellow")

    # silver :
    p7 = p7.difference(p8)
    x, y = p7.exterior.xy  # gold
    plt.fill(x, y, color="gray")

    C = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    p = shapely.Polygon(C).intersection(domain)
    x, y = p.exterior.xy
    plt.fill(x, y, color="brown")

    E = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    p = shapely.Polygon(E)
    x, y = p.exterior.xy
    plt.fill(x, y, color="gray")

    E = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    p = shapely.Polygon(E)
    x, y = p.exterior.xy
    plt.fill(x, y, color="gray")

    E = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    p = shapely.Polygon(E)
    x, y = p.exterior.xy
    plt.fill(x, y, color="gray")

    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    p = shapely.Polygon(E)
    x, y = p.exterior.xy
    plt.fill(x, y, color="gray")

    plt.show()


if __name__ == "__main__":

    # test2()
    # exit()

    import matplotlib.pyplot as plt
    R = Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    C = Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)

    C_cutout = shapely.Polygon(Circle(midpoint=(0. ,0.75), radius = 0.2).discretize_hole(refs=50))

    R_shapely = shapely.Polygon(R)
    C_shapely = shapely.Polygon(C)
    domain = R_shapely.union(C_shapely).difference(C_cutout)
    topo = Topology(domain)

    # P2
    p2 = shapely.Polygon(
        Circle(midpoint=(-0.25, 0.), radius=0.25).discretize_hole(refs=50)
    )
    topo.add_bubble(p2, level=10, is_hole=True)



  
    C_touch_both = shapely.Polygon(Circle(midpoint=(-0.1 ,0.4), radius = 0.3).discretize_hole(refs=50))
    topo.add_bubble(C_touch_both, level=11, is_hole=False)

    # P3
    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    p3 = shapely.Polygon(C3)
    topo.add_bubble(p3, level=1, is_hole=False)

  
    # P4
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    p4 = shapely.Polygon(C4)
    topo.add_bubble(p4, level=1, is_hole=False)

    print(f"Number of bubbles: {len(topo.bubbles)}")


    # this guy is intersection with the domain
    E1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize_hole(refs=50)
    # this guy is just interior
    E2 = Ellipse(midpoint=(-0.55, -0.6), axis=(0.25, 0.1), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E1), level=1, is_hole=False)
    topo.add_bubble(shapely.Polygon(E2), level=1, is_hole=False)




    C5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    C6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    topo.add_bubble(shapely.Polygon(C5), level=3, is_hole=False)
    topo.add_bubble(shapely.Polygon(C6), level=2, is_hole=False)

    C = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    topo.add_bubble(shapely.Polygon(C), level=1, is_hole=False)


 


    E = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)

    E = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)

    E = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)



    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)

    topo.plot_setup()

    quit()

    print(f"Number of bubbles: {len(topo.bubbles)}")

    P = sorted(topo.bubbles, key=lambda p: p.level)
    # check order
    lvl2cl = {0: "blue", 1: "brown", 2: "gray", 3: "yellow", 10: "white"}
    for p in P:
        if isinstance(p.polygon, shapely.MultiPolygon):
            continue
            for p_sub in p.polygon.geoms:
                x, y = p_sub.exterior.xy
        if isinstance(p.polygon, shapely.GeometryCollection):
            continue
        if isinstance(p.polygon, shapely.MultiLineString):
            continue
        print(type(p.polygon))
        x, y = p.polygon.exterior.xy  # gold
        if p.is_hole:
            plt.fill(x, y, color=lvl2cl[p.level])
        else:
            plt.fill(x, y, color=lvl2cl[p.level])

    plt.show()
