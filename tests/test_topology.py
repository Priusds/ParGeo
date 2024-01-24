from bubbles.geometry import Circle, Rectangle
from bubbles.topology import Topology


def test_general():
    # TODO: Add tests
    # Create a topology
    topo = Topology(domain=Rectangle(midpoint=(0,0), width=1, height=1).discretize(), holes={1})

    # Add a circle
    circle = Circle(midpoint=(0.5, 0.5), radius=0.2).discretize(refs=10)
    topo.add(polygon=circle, level=1)