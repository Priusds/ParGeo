import math

from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle, Stellar
from pargeo.gmsh_api import write_geo
from pargeo.transform import Periodic


def generate_domain():
    domain = Circle(midpoint=(0.0, 0.0), radius=1).discretize(
        refs=50
    )  # .difference(Circle(midpoint=(0., 0.), radius = 0.3).discretize(refs=50))
    # domain =  Rectangle(midpoint=(0., 0.), width=1, height=1).to_polygon()
    domain = Domain(domain)

    r11 = Stellar(midpoint=(0.0, 0.0), radius=0.025).discretize(refs=100)
    # Rectangle(midpoint=(0., 0.), width=0.05, height=0.05).to_polygon()
    r12 = Rectangle(midpoint=(0.05, 0.05), width=0.05, height=0.05).to_polygon()

    # r21 = Rectangle(midpoint=(0.05, 0.), width=0.05, height=0.05).to_polygon()
    # r22 = Rectangle(midpoint=(0., 0.05), width=0.05, height=0.05).to_polygon()

    transform = Periodic("any", 0.1, 0.2, alpha=math.pi / 3)
    domain.add_subdomain(r11, level=2, transform=transform)
    domain.add_subdomain(r12, level=1, transform=transform)
    # domain.add(r21, level=2, transform=transform)
    # domain.add(r22, level=2, transform=transform)

    return domain


def generate_football_domain():
    domain = Circle(midpoint=(0.0, 0.0), radius=1).discretize(refs=50)

    domain = Domain(domain)

    c1 = Circle(midpoint=(0.0, 0.0), radius=0.1).discretize(refs=6)
    c2 = Circle(midpoint=(0.125, 0.125), radius=0.1).discretize(refs=6)

    transform = Periodic()
    transform.__set_periodicty("any", 0.25, 0.25)
    domain.add_subdomain(c1, level=2, transform=transform)

    domain.add_subdomain(c2, level=1, transform=transform)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot()
    write_geo(domain=domain, file_name="periodic_examples", correct_curve_loops=True)
