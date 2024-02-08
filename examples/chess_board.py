from bubbles.geometry import Rectangle
from bubbles.gmsh_utils import write_geo
from bubbles.topology import Topology


def generate_topo():
    """Generate a topology with a chess board pattern."""
    N = 3
    R = Rectangle(midpoint=((N - 1) / 2, (N - 1) / 2), width=N, height=N).to_polygon()

    topo = Topology(R, holes={1})
    for i in range(N):
        for j in range(N):
            R = Rectangle(midpoint=(i, j), width=1.0, height=1.0).to_polygon()
            topo.add(polygon=R, level=(i + j) % 2 + 1)

    return topo


if __name__ == "__main__":
    topo = generate_topo()

    topo.plot()

    write_geo(topology=topo, file_name="chess_board", correct_curve_loops=True)
