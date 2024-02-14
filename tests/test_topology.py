"""Unit tests for the domain module."""
from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle


def test_simple():
    domain = Rectangle(midpoint=(0, 0), width=4, height=4).discretize()
    topo = Domain(background=domain)

    c = Circle(midpoint=(0, 0), radius=1).discretize(refs=4)
    r = Rectangle(midpoint=(1.5, 0.0), width=1, height=0.5).discretize()

    topo.add_subdomain(subdomain=c, level=2)
    topo.add_subdomain(subdomain=r, level=1)
