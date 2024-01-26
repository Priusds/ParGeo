from bubbles.geometry import Circle, Rectangle
from bubbles.gmsh_utils import topology_to_gmsh_entities, write_geo
from bubbles.topology import Topology


def generate_topo():
    d_1 = Rectangle(midpoint=(0.0, 0), width=2, height=2).to_polygon()

    dh_1 = Circle(midpoint=(0.0, 0.0), radius=0.5).discretize(refs=50)

    dh_2 = Circle(midpoint=(0.0, 0.0), radius=0.2).discretize(refs=50)

    dh_3 = Circle(midpoint=(0.0, 0.0), radius=0.125).discretize(refs=50)

    cutout = Circle(midpoint=(0.0, 0.0), radius=0.25).discretize(refs=50)

    domain = d_1
    topo = Topology(domain)

    disc_1 = dh_1.difference(cutout)

    disc_2 = dh_2.difference(dh_3)

    topo.add(polygon=disc_2, level=1)
    topo.add(polygon=disc_1, level=2)

    return topo


if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()

    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="edge_cases", correct_curve_loops=True
    )
