from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle


def generate_domain():
    d_1 = Rectangle(midpoint=(0.0, 0), width=1.1, height=1.1).to_polygon()
    domain = d_1
    domain = Domain(domain)

    c1 = Circle(midpoint=(0.0, 0.0), radius=0.5).discretize(refs=50)
    c3 = Circle(midpoint=(0.0, 0.0), radius=0.2).discretize(refs=50)
    c4 = Circle(midpoint=(0.0, 0.0), radius=0.125).discretize(refs=50)
    c2 = Circle(midpoint=(0.0, 0.0), radius=0.25).discretize(refs=50)
    disc_1 = c1 - c2
    disc_2 = c3 - c4

    domain.add_subdomain(subdomain=disc_2, level=1)
    domain.add_subdomain(subdomain=disc_1, level=2)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Disc Edge Case")
