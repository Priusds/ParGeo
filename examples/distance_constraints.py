import random

import shapely

from pargeo.constraint import DistanceConstraint
from pargeo.domain import Domain
from pargeo.geometry import Rectangle, Stellar


def generate_domain():
    """Generate a rectangular domain with many stellar inclusions."""
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()

    P = shapely.Point(0.5, 0.5)

    domain = Domain(domain)
    constraint = DistanceConstraint()
    constraint.set_distance(1, "all", 0.1, to_boundary=True)
    constraint.set_distance("all", P, 0.2, to_boundary=True)

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

    domain.plot("Distance Constraints")
