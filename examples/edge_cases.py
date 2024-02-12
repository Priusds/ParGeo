import math

from pargeo.domain import Domain
from pargeo.geometry import Circle, Ellipse, Rectangle, Stellar
from pargeo.gmsh_utils import write_geo


def generate_domain():
    """Generate a domain with a few edge cases."""
    # Make the domain
    d_1 = Rectangle(midpoint=(0.0, 0), width=2, height=2).to_polygon()

    d_2 = Circle(midpoint=(1.0, 0), radius=0.5).discretize(refs=50)

    dh_1 = Circle(midpoint=(0.0, 0.75), radius=0.2).discretize(refs=50)

    dh_2 = Circle(midpoint=(0.0, 0.0), radius=0.2).discretize(refs=50)

    dh_3 = Ellipse(midpoint=(0.8, 0.0), axis=(0.05, 0.3), angle=0).discretize(refs=50)

    domain = d_1
    domain = domain.union(d_2)
    domain = domain.difference(dh_1)
    domain = domain.difference(dh_2)
    domain = domain.difference(dh_3)
    domain = Domain(domain, holes={3, 10})

    c_1 = Circle(midpoint=(-0.25, 0.0), radius=0.25).discretize(refs=50)
    domain.add(subdomain=c_1, level=10)

    c_2 = Circle(midpoint=(-0.1, 0.4), radius=0.3).discretize(refs=50)
    domain.add(subdomain=c_2, level=11)

    c_3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize(refs=50)
    domain.add(subdomain=c_3, level=1)

    c_4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize(refs=50)
    domain.add(subdomain=c_4, level=1)

    # this guy is intersection with the domain
    e_1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize(refs=200)
    domain.add(subdomain=e_1, level=1)

    s_1 = Stellar(midpoint=(-0.75, -0.6), radius=0.2).discretize(refs=100)
    domain.add(subdomain=s_1, level=1)

    s_2 = Stellar(midpoint=(-0.5, -0.6), radius=0.25).discretize(refs=100)
    domain.add(subdomain=s_2, level=1)

    c_5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize(refs=50)
    domain.add(subdomain=c_5, level=3)

    c_6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize(refs=50)
    domain.add(subdomain=c_6, level=2)

    c_7 = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize(refs=50)
    domain.add(subdomain=c_7, level=1)

    level = 3
    e_2 = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize(
        refs=50
    )
    domain.add(subdomain=e_2, level=level)

    e_3 = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize(refs=50)
    domain.add(subdomain=e_3, level=level)

    e_4 = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize(refs=50)
    domain.add(subdomain=e_4, level=level)

    # this guy encloses the loop !
    e_5 = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize(refs=50)
    domain.add(subdomain=e_5, level=level)

    c_8 = Circle(midpoint=(-0.5, 0.4), radius=0.45).discretize(refs=50)
    domain.add(subdomain=c_8, level=3)

    c_9 = Circle(midpoint=(0.0, 0.75), radius=0.01).discretize(refs=50)
    domain.add(subdomain=c_9, level=3)

    c_10 = Circle(midpoint=(1, -1.0), radius=0.2).discretize(refs=50)
    domain.add(subdomain=c_10, level=10)

    c_11 = Circle(midpoint=(1, -0.8), radius=0.1).discretize(refs=50)
    domain.add(subdomain=c_11, level=2)

    return domain


if __name__ == "__main__":
    domain = generate_domain()

    domain.plot("Edge Cases")

    write_geo(domain=domain, file_name="edge_cases", correct_curve_loops=True)
