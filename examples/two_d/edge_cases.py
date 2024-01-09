import shapely
import math
from bubbles.gmsh_api import write_geo
from bubbles.two_d.hole import Circle, Rectangular, Stellar, Ellipse
from bubbles.two_d.topo_to_gmsh import topology_to_gmsh_entities
from bubbles.two_d.topology_v2 import Topology, Bubble


def generate_topo():
    """Generate a topology with a few edge cases."""
    C = Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)

    C_cutout = shapely.Polygon(
        Circle(midpoint=(0.0, 0.75), radius=0.2).discretize_hole(refs=50)
    )

    R_shapely = shapely.Polygon(
        Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    )
    C_shapely = shapely.Polygon(C)
    C_cutout2 = shapely.Polygon(
        Circle(midpoint=(0.0, 0.0), radius=0.2).discretize_hole(refs=50)
    )
    cutout_3 = shapely.Polygon(
            Ellipse(midpoint=(0.8, 0.0), axis=(.05,.3), angle=0).discretize_hole(refs=50) 
        )
    domain = R_shapely
    domain = domain.union(C_shapely)
    domain = domain.difference(C_cutout)
    domain = domain.difference(C_cutout2)
    domain = domain.difference(cutout_3)
    topo = Topology(domain)


    # P2
    c1 = shapely.Polygon(
        Circle(midpoint=(-0.25, 0.0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add_bubble(Bubble(c1, level=10, is_hole=True))


    c2 = shapely.Polygon(
        Circle(midpoint=(-0.1, 0.4), radius=0.3).discretize_hole(refs=50)
    )
    topo.add_bubble(Bubble(c2, level=11, is_hole=False))

    # P3
    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    p3 = shapely.Polygon(C3)
    topo.add_bubble(Bubble(p3, level=1, is_hole=False))

    # # P4
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    p4 = shapely.Polygon(C4)
    topo.add_bubble(Bubble(p4, level=1, is_hole=False))

    # this guy is intersection with the domain
    E1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize_hole(refs=200)
    topo.add_bubble(Bubble(shapely.Polygon(E1), level=1, is_hole=False))


    S1 = Stellar(midpoint=(-0.75, -0.6), radius=0.2).discretize_hole(refs=100)

    S2 = Stellar(midpoint=(-0.5, -0.6), radius=0.25).discretize_hole(refs=100)

    topo.add_bubble(Bubble(shapely.Polygon(S1), level=1, is_hole=False))
    topo.add_bubble(Bubble(shapely.Polygon(S2), level=2, is_hole=False))


    C5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    C6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    topo.add_bubble(Bubble(shapely.Polygon(C5), level=3, is_hole=False))
    topo.add_bubble(Bubble(shapely.Polygon(C6), level=2, is_hole=False))

    if True:
        C = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
        topo.add_bubble(Bubble(shapely.Polygon(C), level=1, is_hole=False))

        is_hole = False
        level = 3
        E = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
            refs=50
        )
        topo.add_bubble(Bubble(shapely.Polygon(E), level=level, is_hole=is_hole))

        E = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
            refs=50
        )
        topo.add_bubble(Bubble(shapely.Polygon(E), level=level, is_hole=is_hole))

        E = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
            refs=50
        )
        topo.add_bubble(Bubble(shapely.Polygon(E), level=level, is_hole=is_hole))

        # this guy encloses the loop !
        E = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
            refs=50
        )
        topo.add_bubble(Bubble(shapely.Polygon(E), level=level, is_hole=is_hole))

    if True:
        C_touch_hole = shapely.Polygon(
            Circle(midpoint=(-0.5, 0.4), radius=0.45).discretize_hole(refs=50)
        )
        topo.add_bubble(Bubble(C_touch_hole, level=3, is_hole=False))

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()

    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="edge_cases", correct_curve_loops=True
    )
