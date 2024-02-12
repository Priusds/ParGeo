from pargeo.domain import Domain
from pargeo.geometry import Circle
from pargeo.gmsh_utils import write_geo


def generate_domain():
    """Generate a domain with a sequence of centered inclusions."""
    C = Circle(midpoint=(0.0, 0), radius=1.0).discretize(refs=50)

    domain = Domain(C)

    radii = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    level = [i + 1 for i in range(len(radii))]

    for lvl, rad in zip(level, radii):
        C = Circle(midpoint=(0, 0), radius=rad).discretize(refs=50)
        domain.add(subdomain=C, level=lvl)

    return domain


if __name__ == "__main__":
    domain = generate_domain()

    domain.plot("Circle Tunnel")

    write_geo(domain, file_name="circle_tunnel", correct_curve_loops=True)
