"""API for Gmsh.

This module provides an API for Gmsh. It allows the user to write a Gmsh file 
and create a mesh using Gmsh.
"""
from pathlib import Path

from pargeo.domain import Domain
from pargeo.utils.gmsh_utils import (
    domain_to_entities,
    write_geo_from_entities,
    write_msh_from_entities,
)


def write_geo(
    domain: Domain,
    file_name: Path | str,
    correct_curve_loops: bool = False,
) -> None:
    """Write a Gmsh file.

    The file will be saved as a `.geo_unrolled` file that can be opened in Gmsh.

    Args:
        domain (Domain): The domain to be written to a file.
        file_name (Path | str): The file name to save the domain to.
        correct_curve_loops (bool, optional): Whether Gmsh should correct the curve loops.
            Defaults to False.
    """
    gmsh_entities = domain_to_entities(domain)
    write_geo_from_entities(Path(file_name), gmsh_entities, correct_curve_loops)


def write_mesh(
    domain: Domain,
    file_name: Path | str,
    write_geo: bool = True,
    correct_curve_loops: bool = False,
    save_all: bool = False,
) -> None:
    """Create and write mesh with Gmsh.

    Use the Gmsh software to create a mesh and write it to a `.MSH` file.

    Args:
        domain (Domain): The domain to be meshed.
        file_name (Path | str): The file name to save the mesh to. Suffix `.msh` will be added.
        write_geo (bool, optional): Whether to write a `.geo_unrolled` file. Defaults to True.
        correct_curve_loops (bool, optional): Whether to correct the curve loops. Defaults to False.
        save_all (bool, optional): Whether to save all files. Defaults to False. TODO: What does this mean?
    """
    gmsh_entities = domain_to_entities(domain)
    write_msh_from_entities(
        gmsh_entities,
        Path(file_name),
        write_geo,
        correct_curve_loops,
        save_all,
    )
