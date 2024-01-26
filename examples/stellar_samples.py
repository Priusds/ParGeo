import random

from bubbles.geometry import Rectangle, Stellar
from bubbles.gmsh_utils import topology_to_gmsh_entities, write_geo
from bubbles.topology import Topology


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).to_polygon()
    topo = Topology(domain)

    n_stellar = 100
    random.seed(0)
    levels = [random.choice([1, 2, 3, 4, 5]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.5)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        topo.add(polygon=C, level=lvl, extend_domain=True)

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities,
        file_name="stellar_samples",
        correct_curve_loops=True,
    )
