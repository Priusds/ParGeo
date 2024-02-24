from pargeo.domain import Domain
from pargeo.geometry import Rectangle


def generate_domain(n_cells: int = 8) -> Domain:
    """Generate a domain with a chess board pattern."""
    R = Rectangle(
        midpoint=((n_cells - 1) / 2, (n_cells - 1) / 2), width=n_cells, height=n_cells
    ).to_polygon()

    domain = Domain(R, holes={1})
    for i in range(n_cells):
        for j in range(n_cells):
            R = Rectangle(midpoint=(i, j), width=1.0, height=1.0).to_polygon()
            domain.add_subdomain(subdomain=R, level=(i + j) % 2 + 1)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Chess Board")
