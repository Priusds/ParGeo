from pargeo.domain import Domain
from pargeo.geometry import Box


def generate_domain(n_cells: int = 8) -> Domain:
    """Generate a domain with a chess board pattern."""
    board = Box.from_center(
        center=((n_cells - 1) / 2, (n_cells - 1) / 2), width=n_cells, height=n_cells
    ).to_polygon()

    domain = Domain(board, holes={1})
    for i in range(n_cells):
        for j in range(n_cells):
            cell = Box.from_center(center=(i, j), width=1.0, height=1.0).to_polygon()
            domain.add_subdomain(subdomain=cell, level=(i + j) % 2 + 1)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Chess Board")
