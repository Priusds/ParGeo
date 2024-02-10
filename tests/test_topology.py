"""Unit tests for the topology module."""
from pargeo.geometry import Circle, Rectangle
from pargeo.topology import Topology


def test_simple():
    domain = Rectangle(midpoint=(0, 0), width=4, height=4).discretize()
    topo = Topology(domain=domain)

    c = Circle(midpoint=(0, 0), radius=1).discretize(refs=4)
    r = Rectangle(midpoint=(1.5, 0.0), width=1, height=0.5).discretize()

    topo.add(polygon=c, level=2)
    topo.add(polygon=r, level=1)
