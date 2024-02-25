import pytest
from pargeo.domain import Domain
from pargeo.geometry import Box


@pytest.fixture
def domain():
    domain = Box.from_center(center=(0, 0), width=2, height=2).to_polygon()
    return Domain(background=domain)
