import random
from bubbles.gmsh_api import topology_to_gmsh_entities, write_geo
from bubbles.geometry import Rectangle, Stellar, Circle
from bubbles.topology import Topology
from bubbles.transform import Periodic

# (
#     Topology,
#     #clip_x,
#     #clip_y,
#     #polygon_distance_req,
#     #boundary_distance_req,
# )


def generate_topo_simple():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).discretize(refs=4)
    topo = Topology(domain)

    c = Circle(midpoint=(0, 0), radius=0.1).discretize(refs=50)

    rand_midpoint = (random.random(), random.random())
    rand_radius = min(0.09, random.random() * 0.5)
    S = Stellar(midpoint=rand_midpoint, radius=rand_radius).discretize(refs=50)

    transform = Periodic()
    transform.set_periodicty("any", 0.3, 0.5)

    topo.add(c, level=1, transform=transform)

    topo.add(S, level=2, transform=transform)

    return topo


def generate_topo_variation():
    domain = Rectangle(midpoint=(0.5, 0.5), width=1, height=1).discretize(refs=4)
    topo = Topology(domain)

    c = Circle(midpoint=(0, 0), radius=0.02).discretize(refs=50)

    rand_midpoint = (random.random(), random.random())
    rand_radius = min(0.05, random.random() * 0.2)
    S = Stellar(midpoint=rand_midpoint, radius=rand_radius).discretize(refs=50)

    transform = Periodic()
    transform.set_periodicty([1, 4], [0.07, 0.1], [0.2, 0.15])

    topo.add(c, level=4, transform=transform)

    topo.add(S, level=1, transform=transform)

    return topo


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""

    W = 1.0
    H = 1.0

    domain = Rectangle(midpoint=(0.5, 0.5), width=W, height=H).discretize(refs=4)

    topo = Topology(domain)

    n_stellar = 100
    # random.seed(0)
    levels = [random.choice([1, 2, 3, 4, 5]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.045, random.random() * 0.25)) for _ in range(n_stellar)]

    transform = Periodic()
    transform.set_periodicty([1, 2, 3, 4, 5], W, H)

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = Stellar(midpoint=mid, radius=rad).discretize(refs=50)
        topo.add(
            polygon=C,
            level=lvl,
            transform=transform,
        )

    return topo


if __name__ == "__main__":
    topo = generate_topo_variation()  # generate_topo_simple() #generate_topo()
    topo.set_holes({2})
    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities,
        file_name="stellar_samples_periodic_constrains",
        correct_curve_loops=True,
    )
