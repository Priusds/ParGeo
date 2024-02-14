import random

import shapely

from pargeo.constraint import DistanceConstraint
from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle, Stellar
from pargeo.gmsh_api import write_geo


def debug_constraint():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    domain = Domain(domain)

    constraint = DistanceConstraint()
    constraint.set_distance(0.1, 1, "any")

    c1 = Circle(midpoint=(0.4, 0.5), radius=0.1).discretize(refs=50)
    domain.add_subdomain(c1, level=1, constraint=constraint)

    c2 = Circle(midpoint=(0.5, 0.5), radius=0.1).discretize(refs=50)
    domain.add_subdomain(c2, level=2, constraint=constraint)

    return domain


def generate_domain():
    """Generate a rectangular domain with many stellar inclusions."""
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=0.7).to_polygon()
    middle_subdomain = Rectangle(
        midpoint=(0.5, 0.5), width=0.5, height=0.3
    ).to_polygon()

    P = shapely.Point(0.5, 0.5)

    domain = Domain(domain)
    constraint = DistanceConstraint()
    # constraint.set_distance(1, "any", 0.05, to_boundary=True)
    constraint.set_distance("any", middle_subdomain, 0.002, to_boundary=True)
    constraint.show()

    # constraint.set_level_distance(0.05, 1, 2)

    n_stellar = 200
    # random.seed(0)
    levels = [random.choice([1, 2, 3]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.1)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        domain.add_subdomain(
            subdomain=C,
            level=lvl,
            transform=None,  # lambda x: clip_x(x, 0, 1),
            constraint=constraint,
        )

    return domain


if __name__ == "__main__":
    domain = generate_domain()

    domain.plot()

    write_geo(
        domain=domain,
        file_name="stellar_samples_constrains",
        correct_curve_loops=True,
    )
