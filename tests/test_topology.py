from bubbles.geometry import Circle, Rectangle
from bubbles.gmsh_utils import topology_to_gmsh_entities, write_geo
from bubbles.topology import Topology


def test_general():
    # TODO: Add tests
    # Create a topology
    topo = Topology(domain=Rectangle(midpoint=(0,0), width=1, height=1).discretize(), holes={1})

    # Add a circle
    circle = Circle(midpoint=(0.5, 0.5), radius=0.2).discretize(refs=10)
    topo.add(polygon=circle, level=1)

def test_all():
    domain = Rectangle(midpoint=(0,0), width=4, height=4).discretize()
    topo = Topology(domain=domain)

    c = Circle(midpoint=(0, 0), radius=1).discretize(refs=4)
    r = Rectangle(midpoint=(1.5, 0.), width=1, height=0.5).discretize()

    topo.add(polygon=c, level=2) 
    topo.add(polygon=r, level=1)
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="chess_board", correct_curve_loops=True
    )
test_all()