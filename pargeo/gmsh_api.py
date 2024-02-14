from pathlib import Path

from pargeo.domain import Domain
from pargeo.utils.gmsh_utils import (domain_to_entities,
                                     write_geo_from_entities,
                                     write_msh_from_entities)


def write_geo(
    domain: Domain,
    file_name: Path | str,
    correct_curve_loops: bool = False,
) -> None:
    """Convert a domain to GmshEntities."""
    gmsh_entities = domain_to_entities(domain)
    write_geo_from_entities(Path(file_name), gmsh_entities, correct_curve_loops)


def write_mesh(
    domain: Domain,
    file_name: Path | str,
    write_geo: bool = True,
    correct_curve_loops: bool = False,
    save_all: bool = False,
) -> None:
    """Mesh the domain using gmsh."""
    gmsh_entities = domain_to_entities(domain)
    write_msh_from_entities(
        gmsh_entities,
        Path(file_name),
        write_geo,
        correct_curve_loops,
        save_all,
    )
