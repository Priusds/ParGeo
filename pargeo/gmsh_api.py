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

