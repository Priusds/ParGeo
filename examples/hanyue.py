import random

from pargeo.constraint import DistanceConstraint
from pargeo.domain import Domain
from pargeo.geometry import Rectangle, Stellar
from pargeo.transform import Periodic


def generate_domain():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    domain = Domain(domain)

    constraint = DistanceConstraint()
    delta = 0.05
    constraint.set_distance("all", "all", delta)

    transform = Periodic("all", 1.0, 1.0, alpha=0.0)

    n_stellar = 200
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.1)) for _ in range(n_stellar)]

    for mid, rad in zip(midpoints, radii):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        domain.add_subdomain(
            subdomain=C,
            level=1,
            transform=transform,
            constraint=constraint,
        )

    return domain


if __name__ == "__main__":
    domain = generate_domain()

    domain.plot()
