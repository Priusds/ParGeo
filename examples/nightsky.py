import math
import random

from pargeo.constraint import DistanceConstraint
from pargeo.domain import Domain
from pargeo.geometry import NStar, RainDrop, Rectangle


def generate_domain(n_stars: int = 100, n_rain: int = 200) -> Domain:
    """Generate a domain with stars and raindrops."""
    domain = Rectangle(midpoint=(0.0, 0.0), width=2, height=1).to_polygon()
    domain = Domain(domain)

    random.random()

    midpoints = [
        (2 * random.random() - 1, random.random() - 0.5) for _ in range(n_stars)
    ]
    alphas = [2 * math.pi * random.random() for _ in range(n_stars)]

    constraint = DistanceConstraint()
    constraint.set_distance("all", "all", 0.01)

    radii_in = [0.04 * random.random() + 0.005 for _ in range(n_stars)]
    radii_out = [3 * r_in for r_in in radii_in]

    Ns = [random.choice([5]) for _ in range(n_stars)]

    for M, r_in, r_out, N, alph in zip(midpoints, radii_in, radii_out, Ns, alphas):
        star = NStar(M, r_in, r_out, N, alph).to_polygon()
        domain.add_subdomain(star, level=1, constraint=constraint)

    midpoints = [
        (2 * random.random() - 1, random.random() - 0.5) for _ in range(n_rain)
    ]
    A = [0.2 * random.random() + 0.6 for _ in range(n_rain)]
    S = [0.001 * random.random() + 0.005 for _ in range(n_rain)]
    for M, a, scale in zip(midpoints, A, S):
        raindrop = RainDrop(M, a, scale).discretize(128)
        domain.add_subdomain(raindrop, level=2, constraint=constraint)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Night Sky")
