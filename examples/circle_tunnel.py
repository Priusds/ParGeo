from pargeo.domain import Domain
from pargeo.geometry import Circle


def generate_domain(n_circs: int = 10, n_refs: int = 50) -> Domain:
    """Generate a domain with a sequence of centered inclusions."""
    domain = Domain(Circle(midpoint=(0.0, 0), radius=1.0).discretize(refs=50))    
    radii = [1 - r / n_circs for r in range(n_circs)]
    levels = [i + 1 for i in range(n_circs)]

    for lvl, rad in zip(levels, radii):
        domain.add_subdomain(
            subdomain=Circle(midpoint=(0, 0), radius=rad).discretize(refs=n_refs),
            level=lvl
        )

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Circle Tunnel")
