import shapely
import math

from bubbles.two_d.hole import Rectangular, Circle, Ellipse

class Polygon(shapely.Polygon):
    def __init__(self, level, shell=None, holes=None):
        self.level = level

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

    def __init__(self):
        # polygons is a list of non-intersecting polygons (0-measure intersection)
        self.polygons = {}

    
        # HashMap:  level to list(polygons) map 
        self.lvl2polygons = {}   # is final result for gmsh

    # @property
    # def domain(self):
    #     return fun(self)
    
    

    def set_domain(self, domain):
        pass
    
    def __add_discretized_void(self, discretized_hole, level):
        pass


    def __add_discretized_inclusion(self, discretized_hole, level, physical_group): 
        pass

    def add_bubble(self, discretize_hole, level, is_hole): 
        # 1. Add hole/inclusion q
        q = shapely.Polygon(level, discretize_hole)

        I_lower = []
        I_higher = []
        I_same = []
        I_remain = []

        # 2. Find intersections I = p1,...,pn (elements from self.polygons)
        for p in self.polygons: 
            if p.intersects(q):
                if level > p.level : 
                    I_lower.append(p)
                elif level < p.level:
                    I_higher.append(p)
                elif p.level == q : 
                    I_same.append(p)
                else:
                    raise ValueError("Given level is not an integer with order relation.")
            else: 
                I_remain.append(p)
   
        # # update the shape of q constrained by higher level bubbles
        # 3. Find pj1, .., pjk from I such that all have higher level 
        
        # 3.1 If q is completely contained within merge(pj1,...,pjk) -> ignore this hole
        p_all_higher = shapely.union_all(I_higher)
        if p_all_higher.contains(q): 
            # do nothing
            return 
        
        # 3.2 modify q = q.difference(pj1, .., pjk) => q might be a multipolygon
        q = q.difference(p_all_higher)

        # # update the lower level bubbles according to the (remaining) q
        # 4.  Find pj_1, .., pj_l from I such that all have lower level 
        #     and modify pj_k = pj_k.difference(q) for all k = 1 .. l  => pj_k might be a multipolygon
        # 4.1 For each original pj_k that is fully contained in q: remove pj_k  ( maybe pj_k is empty set then ? )
         
        # 4.2 Else if pj_k is a Multipolygon then split it into multiple polygons and delete pj_k
        I_lower_updated = []
        for p in I_lower:    
            if q.contains(p): 
                continue
            else: 
                p_diff = p.difference(q)
                for p_sub in p_diff.geoms:
                    I_lower_updated.append(Polygon(p.level, p_sub.shell, p_sub.holes))

        # # At this point q might be a MultiPolygon
        # 5. Find pi1,...pim, subset of I, s.t. all have the same level as q
        # 5.1 Merge pi1,...pim and q. This is a (possible) multipolygon denoted as P.

        P = shapely.union_all(I_same + [q])
        # 5.2 Delete pi1,...pim and add each subpolygon of P and if P is a multipolygon split into polygons.
        I_same = [Polygon(level, p.shell, p.holes) for p in P.geoms]

        self.polygons = I_higher + I_lower_updated + I_same + I_remain


    def add_void(self, discretize_hole, level):
        # assert hole has discretize method
        pass

    def add_inclusion(self, discretize_hole, level): 
        # assert hole has discretize method
        pass

    def set_physical_group(self, level_to_physical_group_map):
        pass

def test(): 

    lvl2cl = { 0 : 'blue', 1 : 'brown', 2 : 'gray', 3 : 'yellow', 10 : 'white'}


    topo = Topology()

    R = Rectangular(midpoint=(0.,0), width = 2, height = 2).discretize_hole(refs = 4)
    C = Circle(midpoint=(1.,0), radius=0.5).discretize_hole(refs = 50) 
    topo.add_bubble(R, level = 0, is_hole = False)
    topo.add_bubble(C, level = 0, is_hole = False)



    C2 = Circle(midpoint=(-0.25,0), radius=0.25).discretize_hole(refs = 50) 
    topo.add_bubble(C2, level = 10, is_hole = True)


    C3 = Circle(midpoint=(0.75,0), radius=0.25).discretize_hole(refs = 50) 
    C4 = Circle(midpoint=(1.,0), radius=0.25).discretize_hole(refs = 50) 

    topo.add_bubble(C3, level = 1, is_hole = False)
    topo.add_bubble(C4, level = 1, is_hole = False)

    # this guy is intersection with the domain
    E1 = Ellipse(midpoint=(-0.5,0.75), axis=(1,0.1), angle = 0.25 * math.pi).discretize_hole(refs = 50) 
    # this guy is just interior
    E2 = Ellipse(midpoint=(-0.55,-0.6), axis=(0.25,0.1), angle = 0.).discretize_hole(refs = 50) 
    topo.add_bubble(E1, level = 1, is_hole = False)
    topo.add_bubble(E2, level = 1, is_hole = False)


    C5 = Circle(midpoint=(0.5,0.6), radius=0.2).discretize_hole(refs = 50) 
    C6 = Circle(midpoint=(0.6,0.6), radius=0.2).discretize_hole(refs = 50) 
    topo.add_bubble(C6, level = 4, is_hole = False)
    topo.add_bubble(C5, level = 3, is_hole = False)



       
    C = Circle(midpoint = (0.24, -0.8), radius = 0.22).discretize_hole(refs = 50)
    topo.add_bubble(C, level = 1, is_hole = False)

    E = Ellipse(midpoint=(0.25,-0.6), axis=(0.05,0.2), angle = 0.).discretize_hole(refs = 50) 
    topo.add_bubble(E, level = 2, is_hole = False)

    E = Ellipse(midpoint=(0.4,-0.4), axis=(0.2,0.05), angle = 0.).discretize_hole(refs = 50) 
    topo.add_bubble(E, level = 2, is_hole = False)

    E = Ellipse(midpoint=(0.5,-0.6), axis=(0.05,0.2), angle = 0.).discretize_hole(refs = 50) 
    topo.add_bubble(E, level = 2, is_hole = False)

    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4,-0.7), axis=(0.2,0.05), angle = 0.).discretize_hole(refs = 50) 
    topo.add_bubble(E, level = 2, is_hole = False)



    P = sorted(topo.polygons, key = lambda p :p.level)
    # check order 

    for p in P : 
        x,y = p.exterior.xy    # gold
        if p.is_hole : 
            plt.fill(x,y, color = lvl2cl[p.level])
        else: 
            plt.fill(x,y, color =  lvl2cl[p.level])
    
    plt.show()





if __name__ == '__main__':

    # test()

    # exit()

    topo = Topology()

    R = Rectangular(midpoint=(0.,0), width = 2, height = 2).discretize_hole(refs = 4)
    C = Circle(midpoint=(1.,0), radius=0.5).discretize_hole(refs = 50) 

    import matplotlib.pyplot as plt

    # x,y = zip(*R)
    # plt.scatter(x,y,s=2, label = "Rectangular")

    # x,y = zip(*C)
    # plt.scatter(x,y,s=2, color = 'r', label = "Circle")

    # plt.show()

    R_shapely = shapely.Polygon(R)
    C_shapely = shapely.Polygon(C)
    domain = R_shapely.union(C_shapely)

    #x,y = domain.exterior.xy
    #plt.plot(x,y) 
    
    
    topo.set_domain(domain)


    C2 = Circle(midpoint=(-0.25,0), radius=0.25).discretize_hole(refs = 50) 
    p2 = shapely.Polygon(C2)

    # topo.add_discretized_void(C2, level = 10)


    C3 = Circle(midpoint=(0.75,0), radius=0.25).discretize_hole(refs = 50) 
    C4 = Circle(midpoint=(1.,0), radius=0.25).discretize_hole(refs = 50) 
    p3 = shapely.Polygon(C3)
    p4 = shapely.Polygon(C4)


    E1 = Ellipse(midpoint=(-0.5,0.75), axis=(1,0.1), angle = 0.25 * math.pi).discretize_hole(refs = 50) 
    p5 = shapely.Polygon(E1)

    E2 = Ellipse(midpoint=(-0.55,-0.6), axis=(0.25,0.1), angle = 0.).discretize_hole(refs = 50) 
    p6 = shapely.Polygon(E2)



    C5 = Circle(midpoint=(0.5,0.6), radius=0.2).discretize_hole(refs = 50) 
    C6 = Circle(midpoint=(0.6,0.6), radius=0.2).discretize_hole(refs = 50) 
    p7 = shapely.Polygon(C5)
    p8 = shapely.Polygon(C6)
    

    # topo.add_discretized_void(C2, level = 10)


    # topo.__add_discretized_inclusion(C3, level = 2, physical_group="copper")
    # topo.__add_discretized_inclusion(C4, level = 2, physical_group="copper")

    # topo.__add_discretized_inclusion(C6, level = 4, physical_group="gold")
    # topo.__add_discretized_inclusion(C5, level = 3, physical_group="silver")


    


    domain = domain.difference(p2)

    x,y = domain.exterior.xy
    plt.fill(x,y, color = 'blue')

    x,y = domain.interiors[0].xy
    plt.fill(x,y, color = 'white')

    inclusion = p3.union(p4)
    x,y = inclusion.exterior.xy
    plt.fill(x,y, color = 'brown')

    p5 = p5.intersection(domain)
    x,y = p5.exterior.xy
    plt.fill(x,y, color = 'brown')

    x,y = p6.exterior.xy
    plt.fill(x,y, color = 'brown')


    x,y = p8.exterior.xy    # gold
    plt.fill(x,y, color = 'yellow')

    #silver : 
    p7 = p7.difference(p8)
    x,y = p7.exterior.xy    # gold
    plt.fill(x,y, color = 'gray')

    
    C = Circle(midpoint = (0.24, -0.8), radius = 0.22).discretize_hole(refs = 50)
    p = shapely.Polygon(C).intersection(domain)
    x,y = p.exterior.xy    
    plt.fill(x,y, color = 'brown')

    E = Ellipse(midpoint=(0.25,-0.6), axis=(0.05,0.2), angle = 0.).discretize_hole(refs = 50) 
    p = shapely.Polygon(E)
    x,y = p.exterior.xy    
    plt.fill(x,y, color = 'gray')

    E = Ellipse(midpoint=(0.4,-0.4), axis=(0.2,0.05), angle = 0.).discretize_hole(refs = 50) 
    p = shapely.Polygon(E)
    x,y = p.exterior.xy    
    plt.fill(x,y, color = 'gray')

    E = Ellipse(midpoint=(0.5,-0.6), axis=(0.05,0.2), angle = 0.).discretize_hole(refs = 50) 
    p = shapely.Polygon(E)
    x,y = p.exterior.xy    
    plt.fill(x,y, color = 'gray')

    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4,-0.7), axis=(0.2,0.05), angle = 0.).discretize_hole(refs = 50) 
    p = shapely.Polygon(E)
    x,y = p.exterior.xy    
    plt.fill(x,y, color = 'gray')


    plt.show()
