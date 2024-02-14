from pargeo.domain import Domain
from pargeo.geometry import Rectangle
from pargeo.gmsh_api import write_geo


def generate_domain():
    """Generate a domain with a chess board pattern."""
    N = 3
    R = Rectangle(midpoint=((N - 1) / 2, (N - 1) / 2), width=N, height=N).to_polygon()

    domain = Domain(R, holes={1})
    for i in range(N):
        for j in range(N):
            R = Rectangle(midpoint=(i, j), width=1.0, height=1.0).to_polygon()
            domain.add_subdomain(subdomain=R, level=(i + j) % 2 + 1)

    return domain


if __name__ == "__main__":
    domain = generate_domain()

    domain.plot("Chess Board")

    write_geo(domain=domain, file_name="chess_board", correct_curve_loops=True)
