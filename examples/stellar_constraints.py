import random

import shapely

from pargeo.constraint import DistanceConstraint
from pargeo.geometry import Circle, Rectangle, Stellar
from pargeo.gmsh_utils import write_geo
from pargeo.topology import Topology


def debug_constraint():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    topo = Topology(domain)

    constraint = DistanceConstraint()
    constraint.set_distance(0.1, 1, "any")

    c1 = Circle(midpoint=(0.4, 0.5), radius=0.1).discretize(refs=50)
    topo.add(c1, level=1, constraint=constraint)

    c2 = Circle(midpoint=(0.5, 0.5), radius=0.1).discretize(refs=50)
    topo.add(c2, level=2, constraint=constraint)

    return topo


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()

    P = shapely.Point(0.5, 0.5)

    topo = Topology(domain)
    constraint = DistanceConstraint()
    constraint.set_distance(0.05, [1], ["any"], to_boundary=True)
    constraint.set_distance(0.2, "any", P)

    # constraint.set_level_distance(0.05, 1, 2)

    n_stellar = 200
    # random.seed(0)
    levels = [random.choice([1, 2, 3]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.1)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        topo.add(
            polygon=C,
            level=lvl,
            transform=None,  # lambda x: clip_x(x, 0, 1),
            constraint=constraint,
        )

    return topo


if __name__ == "__main__":
    topo = generate_topo()

    topo.plot()

    write_geo(
        topology=topo,
        file_name="stellar_samples_constrains",
        correct_curve_loops=True,
    )
