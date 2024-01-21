import random
import shapely
from bubbles.gmsh_api import topology_to_gmsh_entities, write_geo
from bubbles.two_d.hole import Rectangular, Stellar
from bubbles.two_d.topology import (
    Topology,
    clip_x,
    clip_y,
    polygon_distance_req,
    boundary_distance_req,
)


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""
    domain = shapely.Polygon(
        Rectangular(midpoint=(0.5, 0.5), width=1, height=1).discretize_hole(refs=4)
    )
    topo = Topology(domain)

    n_stellar = 20
    # random.seed(0)
    levels = [random.choice([1, 2, 3, 4, 5]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.5)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = shapely.Polygon(Stellar(midpoint=mid, radius=rad).discretize_hole(refs=50))
        topo.add(
            polygon=C,
            level=lvl,
            transform=lambda x: clip_x(x, 0, 1),
            constraints=(lambda x: boundary_distance_req(x, topo.topology[1], 0.1)) if 1 in topo.topology else None,
        )

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities,
        file_name="stellar_samples_periodic_constrains",
        correct_curve_loops=True,
    )
