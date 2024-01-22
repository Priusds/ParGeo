import shapely
import math
from bubbles.gmsh_api import write_geo, topology_to_gmsh_entities
from bubbles.two_d.hole import Circle, Rectangular, Stellar, Ellipse
from bubbles.two_d.topology_deprecated import Topology_v2 as Topology
from bubbles.two_d.topology_deprecated import Bubble_v2 as Bubble


def generate_topo():
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
    topo = Topology(clip_polygon=domain)
    topo.add(
        polygon=domain, level=0, is_hole=False
    )
    
    # Add bubbles
    c_1 = shapely.Polygon(
        Circle(midpoint=(-0.25, 0.0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_1, level=10, is_hole=True))
    
    
    c_2 = shapely.Polygon(
        Circle(midpoint=(-0.1, 0.4), radius=0.3).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_2, level=11, is_hole=False))
    

    c_3 = shapely.Polygon(
        Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_3, level=1, is_hole=False))
    

    c_4 = shapely.Polygon(
        Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_4, level=1, is_hole=False))
    
    # this guy is intersection with the domain
    e_1 = shapely.Polygon(
        Ellipse(
            midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
        ).discretize_hole(refs=200)
    )
    topo.add(Bubble(polygon=e_1, level=1, is_hole=False))
    

    s_1 = shapely.Polygon(
        Stellar(midpoint=(-0.75, -0.6), radius=0.2).discretize_hole(refs=100)
    )
    topo.add(Bubble(polygon=s_1, level=1, is_hole=False))

    s_2 = shapely.Polygon(
        Stellar(midpoint=(-0.5, -0.6), radius=0.25).discretize_hole(refs=100)
    )
    topo.add(Bubble(polygon=s_2, level=1, is_hole=False))
    

    c_5 = shapely.Polygon(
        Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_5, level=3, is_hole=False))

    c_6 = shapely.Polygon(
        Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_6, level=2, is_hole=False))
    
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

    c_8 = shapely.Polygon(
        Circle(midpoint=(-0.5, 0.4), radius=0.45).discretize_hole(refs=50)
    )
    topo.add(Bubble(polygon=c_8, level=3, is_hole=False))

    c_9 = shapely.Polygon(
        Circle(midpoint=(0.0, 0.75), radius=0.01).discretize_hole(refs=50)
    )

    topo.add(Bubble(polygon=c_9, level=3, is_hole=False))

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    exit()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="edge_cases_v2", correct_curve_loops=True
    )
