from pargeo.domain import Domain
from pargeo.geometry import Circle, Box


def generate_domain():
    """Generate a domain with two discs."""
    d_1 = Box((-1, -1), (1, 1)).to_polygon()
    domain = d_1
    domain = Domain(domain)

    c1 = Circle(center=(0.0, 0.0), radius=0.8).discretize(refs=50)
    c2 = Circle(center=(0.0, 0.0), radius=0.5).discretize(refs=50)
    disc_1 = c1 - c2

    c3 = Circle(center=(0.0, 0.0), radius=0.4).discretize(refs=50)
    c4 = Circle(center=(0.0, 0.0), radius=0.2).discretize(refs=50)
    disc_2 = c3 - c4

    domain.add_subdomain(subdomain=disc_2, level=1)
    domain.add_subdomain(subdomain=disc_1, level=2)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Disc Edge Case")
