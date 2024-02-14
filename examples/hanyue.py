import random

from pargeo.constraint import DistanceConstraint
from pargeo.domain import Domain
from pargeo.geometry import Rectangle, Stellar
from pargeo.gmsh_api import write_geo
from pargeo.transform import Periodic


def generate_domain():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    domain = Domain(domain)

    constraint = DistanceConstraint()
    delta = 0.05
    constraint.set_distance("any", "any", delta)

    # Periodic repeat of geometries of any level with periodic length 1 in both directions
    # here x,y axis direction since alpha = 0.
    transform = Periodic("any", 1.0, 1.0, alpha=0.0)

    # f( x) = f( x + L )

    n_stellar = 1000
    # random.seed(0)
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


def subdomain():
    domain = Rectangle(midpoint=(0.25, 0.25), width=0.5, height=0.5).to_polygon()
    domain = Domain(domain)

    domain_whole = generate_domain()

    lvl2polygon = domain_whole.level_to_subdomain

    level = 1
    domain.add_subdomain(lvl2polygon[level], level)

    import matplotlib.pyplot as plt

    plt.subplot(1, 2, 1)
    domain_whole.plot()

    plt.plot([0.5, 0.5], [0, 1], "-k", linewidth=2)
    plt.plot([0, 1], [0.5, 0.5], "-k", linewidth=2)

    plt.subplot(1, 2, 2)
    domain.plot()
    plt.show()


if __name__ == "__main__":
    # subdomain()
    # exit()
    domain = generate_domain()

    domain.plot()

    write_geo(
        domain=domain,
        file_name="hanyue",
        correct_curve_loops=True,
    )
