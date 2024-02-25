from pargeo.constraint import DistanceConstraint
from pargeo.geometry import Box


def test_distance_constraint(domain):
    constraint = DistanceConstraint()
    constraint.set_distance(obj_1=1, obj_2=1, distance=0.1)
    constraint.set_distance(obj_1=3, obj_2="all", distance=0.1)

    subdomain1 = Box.from_center(center=(0, 0), width=0.1, height=0.1).to_polygon()
    subdomain2 = Box.from_center(
        center=(0.05, 0.05), width=0.1, height=0.1
    ).to_polygon()
    subdomain3 = Box.from_center(
        center=(0.05, 0.05), width=0.1, height=0.1
    ).to_polygon()
    subdomain4 = Box.from_center(
        center=(0.06, 0.04), width=0.1, height=0.1
    ).to_polygon()

    domain.add_subdomain(subdomain1, 1)
    assert not domain.add_subdomain(
        subdomain2, 1, constraint=constraint
    ), "Subdomain 2 should not be added."
    assert len(domain.subdomains) == 2

    assert not domain.add_subdomain(
        subdomain4, 3, constraint=constraint
    ), "Subdomain 4 should not be added."
    assert len(domain.subdomains) == 2

    assert domain.add_subdomain(
        subdomain3, 2, constraint=constraint
    ), "Subdomain 3 should now be added."
    assert len(domain.subdomains) == 3

    assert not domain.add_subdomain(
        subdomain4, 3, constraint=constraint
    ), "Subdomain 4 should not be added."
    assert len(domain.subdomains) == 3
