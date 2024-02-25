from pargeo.geometry import Box


def test_set_holes(domain):
    new_holes = {1, 2, 3}
    domain.set_holes(new_holes)
    assert domain.holes == new_holes


def test_add_subdomain(domain):
    subdomains = [
        Box.from_center(center=(0.5, 0.5), width=0.1, height=0.1).to_polygon(),
        Box.from_center(center=(-0.5, 0.5), width=0.1, height=0.1).to_polygon(),
        Box.from_center(center=(0.5, -0.5), width=0.1, height=0.1).to_polygon(),
        Box.from_center(center=(-0.5, -0.5), width=0.1, height=0.1).to_polygon(),
    ]
    for i, subdomain in enumerate(subdomains):
        assert domain.add_subdomain(subdomain, i + 1)
    assert len(domain.subdomains) == 5

    subdomain = Box.from_center(center=(0.6, 0.6), width=0.12, height=0.12).to_polygon()
    domain.add_subdomain(subdomain, 1)
    assert len(domain.subdomains) == 5

    subdomain = Box.from_center(
        center=(-0.6, -0.6), width=0.12, height=0.12
    ).to_polygon()
    domain.add_subdomain(subdomain, 1)
    assert len(domain.subdomains) == 7 # 6
