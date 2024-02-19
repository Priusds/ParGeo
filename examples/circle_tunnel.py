from pargeo.domain import Domain
from pargeo.geometry import Circle


def generate_domain():
    """Generate a domain with a sequence of centered inclusions."""
    C = Circle(midpoint=(0.0, 0), radius=1.0).discretize(refs=50)

    domain = Domain(C)

    n_radii = 10
    diff = 1.0 / (n_radii)
    radii = [1 - r * diff for r in range(n_radii)]
    level = [i + 1 for i in range(n_radii)]

    for lvl, rad in zip(level, radii):
        C = Circle(midpoint=(0, 0), radius=rad).discretize(refs=50)
        domain.add_subdomain(subdomain=C, level=lvl)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Circle Tunnel")
