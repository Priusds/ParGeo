from pargeo.geometry import Circle
from pargeo.gmsh_utils import write_geo
from pargeo.topology import Topology


def generate_topo():
    """Generate a topology with a sequence of centered inclusions."""
    C = Circle(midpoint=(0.0, 0), radius=1.0).discretize(refs=50)

    topo = Topology(C)

    radii = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    level = [i + 1 for i in range(len(radii))]

    for lvl, rad in zip(level, radii):
        C = Circle(midpoint=(0, 0), radius=rad).discretize(refs=50)
        topo.add(polygon=C, level=lvl)

    return topo


if __name__ == "__main__":
    topo = generate_topo()

    topo.plot()

    write_geo(topo, file_name="circle_tunnel", correct_curve_loops=True)
