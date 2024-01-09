import shapely
from bubbles.gmsh_api import write_geo
from bubbles.two_d.hole import Rectangular
from bubbles.two_d.topo_to_gmsh import topology_to_gmsh_entities
from bubbles.two_d.topology_v2 import Topology


def generate_topo():
    N = 10
    R = shapely.Polygon(
        Rectangular(
            midpoint=((N - 1) / 2, (N - 1) / 2), width=N, height=N
        ).discretize_hole(refs=4)
    )
    topo = Topology(R)

    # TODO: this is buggy as above hole not added correctly

    for i in range(N):
        for j in range(N):
            R = shapely.Polygon(
                Rectangular(midpoint=(i, j), width=1.0, height=1.0).discretize_hole(
                    refs=4
                )
            )
            topo.add_bubble(R, level=(i + j) % 2 + 1, is_hole=(i + j) % 2)

    return topo


topo = generate_topo()
topo.plot()

gmsh_entities = topology_to_gmsh_entities(topo)
write_geo(
    gmsh_entities=gmsh_entities, file_name="chess_board", correct_curve_loops=True
)