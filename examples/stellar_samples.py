import random

from bubbles.geometry import Rectangle, Stellar
from bubbles.gmsh_utils import write_geo, mesh
from bubbles.topology import Topology


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""
    domain = Rectangle(midpoint=(0.5, 0.5), width=2, height=1).to_polygon()
    topo = Topology(domain)

    n_stellar = 800
    # random.seed(0)
    levels = [random.choice([1, 2, 3, 4, 5]) for _ in range(n_stellar)]
    midpoints = [(2*random.random()*random.choice([1,-1]), random.random(),) for _ in range(n_stellar)]

    radii = [min((0.09, random.random() * 0.5)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        topo.add(polygon=C, level=lvl, extend_domain=False)

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    write_geo(
        topology=topo,
        file_name="stellar_samples",
        correct_curve_loops=True,
    )
