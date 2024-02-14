from pargeo.domain import Domain
from pargeo.geometry import Circle, Rectangle
from pargeo.gmsh_api import write_geo


def generate_domain():
    d_1 = Rectangle(midpoint=(0.0, 0), width=2, height=2).to_polygon()

    dh_1 = Circle(midpoint=(0.0, 0.0), radius=0.5).discretize(refs=50)

    dh_2 = Circle(midpoint=(0.0, 0.0), radius=0.2).discretize(refs=50)

    dh_3 = Circle(midpoint=(0.0, 0.0), radius=0.125).discretize(refs=50)

    cutout = Circle(midpoint=(0.0, 0.0), radius=0.25).discretize(refs=50)

    domain = d_1
    domain = Domain(domain)

    disc_1 = dh_1.difference(cutout)

    disc_2 = dh_2.difference(dh_3)

    domain.add_subdomain(subdomain=disc_2, level=1)
    domain.add_subdomain(subdomain=disc_1, level=2)

    return domain


if __name__ == "__main__":
    domain = generate_domain()
    domain.plot("Disc Edge Case")

    write_geo(domain=domain, file_name="edge_cases", correct_curve_loops=True)
