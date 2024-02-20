import pytest
from pargeo.domain import Domain
from pargeo.geometry import Rectangle


@pytest.fixture
def domain():
    domain = Rectangle(midpoint=(0, 0), width=2, height=2).to_polygon()
    return Domain(background=domain)
