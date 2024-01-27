import random

from bubbles.constraint import DistanceConstraint
from bubbles.geometry import Circle, Rectangle, Stellar
from bubbles.gmsh_utils import write_geo
from bubbles.topology import Topology
from bubbles.transform import Periodic


def generate_topo():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    topo = Topology(domain)

    constraint = DistanceConstraint()
    delta = 0.05
    constraint.set_distance(delta, "any", "any")

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
        topo.add(
            polygon=C,
            level=1,
            transform=transform,
            constraint=constraint,
        )

    return topo


def subdomain():
    domain = Rectangle(midpoint=(0.25, 0.25), width=0.5, height=0.5).to_polygon()
    topo = Topology(domain)

    topo_whole = generate_topo()

    lvl2polygon = topo_whole.lvl2multipoly

    level = 1
    topo.add(lvl2polygon[level], level)

    import matplotlib.pyplot as plt

    plt.subplot(1, 2, 1)
    topo_whole.plot()

    plt.plot([0.5, 0.5], [0, 1], "-k", linewidth=2)
    plt.plot([0, 1], [0.5, 0.5], "-k", linewidth=2)

    plt.subplot(1, 2, 2)
    topo.plot()
    plt.show()


if __name__ == "__main__":
    # subdomain()
    # exit()
    topo = generate_topo()

    topo.plot()

    write_geo(
        topology=topo,
        file_name="hanyue",
        correct_curve_loops=True,
    )
