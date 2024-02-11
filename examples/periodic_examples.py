import math

from pargeo.geometry import Circle, Rectangle, Stellar
from pargeo.gmsh_utils import write_geo
from pargeo.topology import Topology
from pargeo.transform import Periodic


def generate_topo():
    domain = Circle(midpoint=(0.0, 0.0), radius=1).discretize(
        refs=50
    )  # .difference(Circle(midpoint=(0., 0.), radius = 0.3).discretize(refs=50))
    # domain =  Rectangle(midpoint=(0., 0.), width=1, height=1).to_polygon()
    topo = Topology(domain)

    r11 = Stellar(midpoint=(0.0, 0.0), radius=0.025).discretize(refs=100)
    # Rectangle(midpoint=(0., 0.), width=0.05, height=0.05).to_polygon()
    r12 = Rectangle(midpoint=(0.05, 0.05), width=0.05, height=0.05).to_polygon()

    # r21 = Rectangle(midpoint=(0.05, 0.), width=0.05, height=0.05).to_polygon()
    # r22 = Rectangle(midpoint=(0., 0.05), width=0.05, height=0.05).to_polygon()

    transform = Periodic("any", 0.1, 0.2, alpha=math.pi / 3)
    topo.add(r11, level=2, transform=transform)
    topo.add(r12, level=1, transform=transform)
    # topo.add(r21, level=2, transform=transform)
    # topo.add(r22, level=2, transform=transform)

    return topo


def generate_topo_football():
    domain = Circle(midpoint=(0.0, 0.0), radius=1).discretize(refs=50)

    topo = Topology(domain)

    c1 = Circle(midpoint=(0.0, 0.0), radius=0.1).discretize(refs=6)
    c2 = Circle(midpoint=(0.125, 0.125), radius=0.1).discretize(refs=6)

    transform = Periodic()
    transform.__set_periodicty("any", 0.25, 0.25)
    topo.add(c1, level=2, transform=transform)

    topo.add(c2, level=1, transform=transform)

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    write_geo(topology=topo, file_name="periodic_examples", correct_curve_loops=True)
