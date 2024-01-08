import shapely
import math
import matplotlib.pyplot as plt
from bubbles.two_d.hole import Rectangular, Circle, Ellipse

from collections import namedtuple
from copy import deepcopy

Bubble = namedtuple("Bubble", ["polygon", "level", "is_hole"])

GRID_SIZE = 1e-15

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

    def __init__(self, polygon: shapely.Polygon, grid_size=GRID_SIZE):
        # polygons is a list of non-intersecting polygons (0-measure intersection)
        #polygon = shapely.set_precision(polygon, grid_size)
        self.bubbles = {Bubble(polygon, level=0, is_hole=False)}

        self.ref_domain = Bubble(polygon, level=0, is_hole=False)
        self.grid_size = grid_size
        # TODO: Check if it is enough just to set the grid_size level for all new geometries,
        # and then leave the parameter out for the rest of the code.

        # HashMap:  level to list(polygons) map
        self.lvl2polygons = {}  # is final result for gmsh

    @property
    def domain(self):
        domain_polygons = [p.polygon for p in self.bubbles if p.level == 0]
        return shapely.MultiPolygon(domain_polygons)

    def add_bubble(self, Q: shapely.Polygon, level: int, is_hole: bool):
        if not level > 0:
            raise ValueError("Level must be positive integer.")
        # clip q with respect to the reference domain
        new_bubble = Q.intersection(self.ref_domain.polygon, grid_size=self.grid_size)
        if not isinstance(new_bubble, shapely.Polygon):
            raise ValueError(f"new_bubble is not a Polygon, found type {type(new_bubble)}.")

        # Find intersecting polygons and group by level relation
        lower_bubbles = set()
        higher_bubbles = set()
        equal_bubbles = set()
        for bubble in self.bubbles:
            if bubble.polygon.intersects(new_bubble):
                if level > bubble.level:
                    lower_bubbles.add(bubble)
                elif level < bubble.level:
                    higher_bubbles.add(bubble)
                elif level == bubble.level:
                    equal_bubbles.add(bubble)
                    if not bubble.is_hole == is_hole:
                        raise ValueError(
                            "Bubbles with same level must have same is_hole value."
                        )
                else:
                    raise ValueError(
                        "Given level is not an integer with order relation."
                    )
        if len(higher_bubbles) > 0:
            # Update the shape of new_bubble constrained by higher level bubbles
            # If new_bubble is completely contained within the union of higher level bubbles -> ignore q
            higher_polygons_union = shapely.union_all(
                [b.polygon for b in higher_bubbles]
            )
            if higher_polygons_union.contains(new_bubble):
                return
            new_bubble = new_bubble.difference(higher_polygons_union, grid_size=self.grid_size)
            # new_bubble might be a multipolygon
            if not (isinstance(new_bubble, shapely.Polygon) or isinstance(new_bubble, shapely.MultiPolygon)):
                raise ValueError(f"new_bubble is of type {type(new_bubble)}")

        for lower_bubble in lower_bubbles:
            p_diff = lower_bubble.polygon.difference(new_bubble, grid_size=self.grid_size)
            if not (isinstance(p_diff, shapely.Polygon) or isinstance(p_diff, shapely.MultiPolygon)):
                raise ValueError(f"p_diff is of type {type(p_diff)}")

            if p_diff.area > 0:
                if isinstance(p_diff, shapely.MultiPolygon):
                    for p_sub in p_diff.geoms:
                        assert isinstance(p_sub, shapely.Polygon)
                        self.bubbles.add(
                            Bubble(p_sub, level=lower_bubble.level, is_hole=lower_bubble.is_hole)
                        )
                        # Add intersection points to new_bubble
                        new_bubble = new_bubble.intersection(p_sub.union(new_bubble, grid_size=self.grid_size), grid_size=self.grid_size)
                else:
                    # Add intersection points to new_bubble
                    new_bubble = new_bubble.intersection(p_diff.union(new_bubble, grid_size=self.grid_size), grid_size=self.grid_size)
                    self.bubbles.add(Bubble(p_diff, level=lower_bubble.level, is_hole=lower_bubble.is_hole))
            else:
                # lower_bubble is covered by new_bubble
                if lower_bubble.level == 0:
                    # new_bubble completely contains the domain, this is not allowed
                    return
                else:
                    # Remove lower_bubble from self.bubbles
                    pass
            self.bubbles.remove(lower_bubble)

        if not (isinstance(new_bubble, shapely.Polygon) or isinstance(new_bubble, shapely.MultiPolygon)):
            raise ValueError(f"new_bubble is of type {type(new_bubble)}")
        # Merge q_diff with equal level bubbles, the result still might be a multipolygon
        for lower_bubble in equal_bubbles:
            self.bubbles.remove(lower_bubble)

        q_final = shapely.union_all([new_bubble] + [p.polygon for p in equal_bubbles],grid_size=self.grid_size)
        if not (isinstance(q_final, shapely.Polygon) or isinstance(q_final, shapely.MultiPolygon)):
            raise ValueError(f"q_final is of type {type(q_final)}")
        
        if isinstance(q_final, shapely.MultiPolygon):
            for q_sub in q_final.geoms:
                    self.bubbles.add(
                        Bubble(q_sub, level=level, is_hole=is_hole)
                    )
        else:
            self.bubbles.add(Bubble(q_final, level=level, is_hole=is_hole))


    def plot_setup(self): 
        boundary_color = 'black'
        reference_domain_hole_color = 'white'
        # @pre lvl2cl does not contain boundary_color
        # @pre lvl2cl does not contain reference_domain_hole_color
        lvl2cl = {0: "blue", 1: "brown", 2: "gray", 3: "yellow", 10: "lightgray", 11: "purple"}

        
        setted_labels = []

        for bubble in self.bubbles: 
            p, lvl, is_hole = bubble 

            if len(p.interiors) > 0 : 
                x,y = p.exterior.xy
                plt.fill(x, y, color=lvl2cl[lvl]) 
                if is_hole : 
                    if not boundary_color in setted_labels: 
                        plt.plot(x, y, color=boundary_color, label = 'domain boundary') 
                        setted_labels.append(boundary_color)
                    else: 
                        plt.plot(x, y, color=boundary_color)
                for pp in p.interiors : 
                    x, y = pp.xy  # gol
                    plt.fill(x, y, color='white') 
                    if is_hole : 
                        if not boundary_color in setted_labels: 
                            plt.plot(x, y, color=boundary_color, label = 'domain boundary') 
                            setted_labels.append(boundary_color)
                        else: 
                            plt.plot(x, y, color=boundary_color)

        for bubble in self.bubbles: 
            p, lvl, is_hole = bubble 

            if len(p.interiors) == 0 : 
                x,y = p.exterior.xy
                color = lvl2cl[lvl]
                if not color in setted_labels: 
                    label = f'Incl. with lvl = {lvl}' if not is_hole else f'Hole with lvl = {lvl}'
                    plt.fill(x, y, color=color, label = label)
                    setted_labels.append(color)
                else:
                    plt.fill(x, y, color=color)
                if is_hole : 
                    if not boundary_color in setted_labels: 
                        plt.plot(x, y, color='black', label = 'domain boundary') 
                        setted_labels.append(boundary_color)
                    else:
                        plt.plot(x, y, color=boundary_color)

        
        x,y = self.ref_domain.polygon.exterior.xy
        plt.plot(x, y, color='black')
        for pp in self.ref_domain.polygon.interiors : 
            x, y = pp.xy  # gol

            if not reference_domain_hole_color in setted_labels:
                plt.fill(x, y, color=reference_domain_hole_color, label = 'domain hole') 
                setted_labels.append(reference_domain_hole_color)
            else:
                plt.fill(x, y, color=reference_domain_hole_color)


            if not boundary_color in setted_labels: 
                plt.plot(x, y, color=boundary_color, label = 'domain boundary')
                setted_labels.append(boundary_color)
            else:
                plt.plot(x, y, color=boundary_color)

        plt.legend()
        plt.show()


def test1():
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

    C_touch_hole = shapely.Polygon(Circle(midpoint=(-0.5 ,0.4), radius = 0.45).discretize_hole(refs=50))
    topo.add_bubble(C_touch_hole, level=3, is_hole=False)

    # # P3
    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    p3 = shapely.Polygon(C3)
    topo.add_bubble(p3, level=1, is_hole=False)

  
    # # P4
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    p4 = shapely.Polygon(C4)
    topo.add_bubble(p4, level=1, is_hole=False)

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
    pass
    # test1()
    test2()
