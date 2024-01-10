from bubbles.gmsh_api import topology_to_gmsh_entities, write_geo
from bubbles.two_d.hole import Rectangular, Circle
from bubbles.two_d.topology import PeriodicTopology, Bubble
import shapely

def generate_topo():
    """Generate a simple periodic topology."""
    topo = PeriodicTopology(
        midpoint=(0,0),
        width=1,
        height=1,
        periodic_x=True,
        periodic_y=True,
    )
    c1 = Bubble(
        polygon=shapely.Polygon(Circle(midpoint=(0.5, 0.5), radius=0.1).discretize_hole(refs=50)),
        level=1, 
        is_hole=False
    )
    topo.add(c1)

    c2 = Bubble(
        polygon=shapely.Polygon(Circle(midpoint=(0.5, 0.), radius=0.1).discretize_hole(refs=50)),
        level=1,
        is_hole=False
    )
    topo.add(c2)

    c3 = Bubble(
        polygon=shapely.Polygon(Circle(midpoint=(0., 0.7), radius=0.1).discretize_hole(refs=50)),
        level=1, 
        is_hole=False
    )
    topo.add(c3)
    return topo

if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="simple_periodicity", correct_curve_loops=True
    )