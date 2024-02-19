import math
import random

from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle, Stellar
from pargeo.gmsh_api import write_geo
from pargeo.transform import Periodic


def generate_simple_domain():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    domain = Domain(domain)

    c = Circle(midpoint=(0, 0), radius=0.1).discretize(refs=50)

    rand_midpoint = (random.random(), random.random())
    rand_radius = min(0.09, random.random() * 0.5)
    S = Stellar(midpoint=rand_midpoint, radius=rand_radius).discretize(refs=50)

    transform = Periodic("all", 0.3, 0.5)

    domain.add_subdomain(c, level=1, transform=transform)

    domain.add_subdomain(S, level=2, transform=transform)

    return domain


def generate_variation_domain():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    domain = Domain(domain)

    c = Circle(midpoint=(0, 0), radius=0.02).discretize(refs=50)

    rand_midpoint = (random.random(), random.random())
    rand_radius = min(0.05, random.random() * 0.2)
    S = Stellar(midpoint=rand_midpoint, radius=rand_radius).discretize(refs=50)

    transform = Periodic([1, 4], [0.07, 0.1], [0.2, 0.15], alpha=math.pi / 4)

    domain.add_subdomain(c, level=4, transform=transform)

    domain.add_subdomain(S, level=1, transform=transform)

    return domain


def generate_domain():
    """Generate a rectangular domain with many stellar inclusions."""

    W = 1.0
    H = 1.0

    domain = Rectangle(midpoint=(0.5, 0.5), width=W, height=H).to_polygon()

    domain = Domain(domain)

    n_stellar = 100
    # random.seed(0)
    levels = [random.choice([1, 2, 3, 4, 5]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.045, random.random() * 0.25)) for _ in range(n_stellar)]

    transform = Periodic()
    transform.__set_periodicty([1, 2, 3, 4, 5], W, H, alpha=math.pi / 4)

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        domain.add_subdomain(
            subdomain=C,
            level=lvl,
            transform=transform,
        )

    return domain


if __name__ == "__main__":
    domain = generate_variation_domain()  # generate_simple_domain() #generate_domain()
    domain.set_holes({2})
    domain.plot()

    write_geo(
        domain=domain,
        file_name="stellar_samples_periodic_constrains",
        correct_curve_loops=True,
    )
