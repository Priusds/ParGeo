import shapely
import math
from bubbles.gmsh_api import bubbles_to_gmsh_entities, write_geo, topology_to_gmsh_entities
from bubbles.two_d.hole import Circle, Rectangular, Stellar, Ellipse
from bubbles.two_d.topology import Topology_v3 as Topology
from bubbles.two_d.topology import Bubble


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
    domain = domain.difference(dh_2)
    #domain = domain.difference(dh_3)
    topo = Topology(mask=domain)
#    topo.add(Bubble(polygon=domain, level=0, is_hole=False))
    topo.add(multi_polygon=domain, level=0)
    
    # Add bubbles
    c_1 = shapely.Polygon(
        Circle(midpoint=(-0.25, 0.0), radius=0.25).discretize_hole(refs=50)
    )
    #topo.add(Bubble(polygon=c_1, level=10, is_hole=True))
    topo.add(multi_polygon=c_1, level=10)
    return topo
    c_2 = shapely.Polygon(
        Circle(midpoint=(-0.1, 0.4), radius=0.3).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_2, level=11)
    

    c_3 = shapely.Polygon(
        Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_3, level=1)

    c_4 = shapely.Polygon(
        Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_4, level=1)
    
    # this guy is intersection with the domain
    e_1 = shapely.Polygon(
        Ellipse(
            midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
        ).discretize_hole(refs=200)
    )
    topo.add(multi_polygon=e_1, level=1)

    s_1 = shapely.Polygon(
        Stellar(midpoint=(-0.75, -0.6), radius=0.2).discretize_hole(refs=100)
    )
    topo.add(multi_polygon=s_1, level=1)

    s_2 = shapely.Polygon(
        Stellar(midpoint=(-0.5, -0.6), radius=0.25).discretize_hole(refs=100)
    )
    topo.add(multi_polygon=s_2, level=1)
    

    c_5 = shapely.Polygon(
        Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_5, level=3)

    c_6 = shapely.Polygon(
        Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_6, level=2)

    c_7 = shapely.Polygon(
        Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    )
    topo.add(multi_polygon=c_7, level=1)
    return topo
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
    from bubbles.two_d.topology import plot2
    topo = generate_topo()

    cut_out_levels = {}

    level2is_hole = dict()
    for lvl in topo.topology.keys():
        level2is_hole[lvl] = True if lvl in cut_out_levels else False
    bubbles = topo.get_bubbles(cut_out_levels=cut_out_levels)

    plot2(bubbles, level2is_hole)

    gmsh_entities = bubbles_to_gmsh_entities(bubbles)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="edge_cases_2", correct_curve_loops=True
    )
