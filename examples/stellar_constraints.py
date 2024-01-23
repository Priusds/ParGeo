import random
import shapely
from bubbles.gmsh_api import topology_to_gmsh_entities, write_geo
from bubbles.two_d.hole import Rectangular, Stellar, Circle
from bubbles.two_d.constraints import DistanceConstraint


from bubbles.two_d.topology import Topology


def debug_constraint():
    domain = shapely.Polygon(
        Rectangular(midpoint=(0.5, 0.5), width=1, height=1).discretize_hole(refs=4)
    )

    topo = Topology(domain)

    c1 = shapely.Polygon(
        Circle(midpoint=(0.4, 0.5), radius=0.1).discretize_hole(refs=50)
    )
    c2 = shapely.Polygon(
        Circle(midpoint=(0.5, 0.5), radius=0.1).discretize_hole(refs=50)
    )
    constraint = DistanceConstraint()
    constraint.set_level_distance(0.1, 1, "any")

    topo.add(c1, level=1, constraints=constraint)

    topo.add(c2, level=2, constraints=constraint)

    return topo


def generate_topo():
    """Generate a rectangular topology with many stellar inclusions."""
    domain = shapely.Polygon(
        Rectangular(midpoint=(0.5, 0.5), width=1, height=1).discretize_hole(refs=4)
    )

    # c1 = shapely.Polygon(Circle(midpoint=(0., 0.), radius=0.1).discretize_hole(refs=50))
    P = shapely.Point(0.5, 0.5)

    topo = Topology(domain)
    constraint = DistanceConstraint()
    constraint.set_level_distance(0.02, 1, "any")
    constraint.set_geometry_distance(P, 0.2, "any")
    # constraint.set_level_distance(0.05, 1, 2)

    n_stellar = 2000
    # random.seed(0)
    levels = [random.choice([1, 2, 3]) for _ in range(n_stellar)]
    midpoints = [(random.random(), random.random()) for _ in range(n_stellar)]
    radii = [min((0.09, random.random() * 0.1)) for _ in range(n_stellar)]

    for mid, rad, lvl in zip(midpoints, radii, levels):
        C = shapely.Polygon(Stellar(midpoint=mid, radius=rad).discretize_hole(refs=50))
        topo.add(
            polygon=C,
            level=lvl,
            transform=None,  # lambda x: clip_x(x, 0, 1),
            constraints=constraint,
        )

    return topo


if __name__ == "__main__":
    topo = generate_topo()

    topo.plot()
    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities,
        file_name="stellar_samples_constrains",
        correct_curve_loops=True,
    )
