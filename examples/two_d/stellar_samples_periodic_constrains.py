import random
import shapely
from bubbles.gmsh_api import topology_to_gmsh_entities, write_geo
from bubbles.two_d.hole import Rectangular, Stellar, Circle
from bubbles.two_d.topology import (
    Topology,
    clip_x,
    clip_y,
    polygon_distance_req,
    boundary_distance_req,
)

from bubbles.two_d.transform import Periodic_New

def generate_topo_simple(): 
    domain = shapely.Polygon(
        Rectangular(midpoint=(0.5, 0.5), width=1, height=1).discretize_hole(refs=4)
    )
    topo = Topology(domain)

    c = shapely.Polygon(
        Circle(midpoint=(0, 0), radius=0.1).discretize_hole(refs=50)
    )

    transform = Periodic_New()
    transform.set_periodicty(1, 0.3, 0.5) 

    topo.add(c, level = 1, transform = transform)

    return topo


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
            constraints=(lambda x: polygon_distance_req(x, topo.topology[1], 0.1)) if 1 in topo.topology else None,
        )

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.set_holes({2,3})
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities,
        file_name="stellar_samples_periodic_constrains",
        correct_curve_loops=True,
    )
