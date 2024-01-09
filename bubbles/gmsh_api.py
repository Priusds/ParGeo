"""Module for mesh creation using gmsh."""
from enum import IntEnum
from pathlib import Path
from typing import Optional, Sequence

import gmsh
from pydantic import BaseModel, StrictInt


class PhysicalDimension(IntEnum):
    """Dimensions allowed for physical groups and/or meshes."""

    two = 2
    three = 3


class Point(BaseModel):
    """Two dimensional point."""

    x: float
    y: float
    z: float
    lc: float
    tag: StrictInt

    class Config:
        allow_mutation = False


class Line(BaseModel):
    """Line between two points."""

    start_tag: StrictInt
    end_tag: StrictInt
    tag: StrictInt

    class Config:
        allow_mutation = False


class CurveLoop(BaseModel):
    """Closed loop of lines."""

    line_tags: Sequence[StrictInt]
    tag: StrictInt

    class Config:
        allow_mutation = False


class PlaneSurface(BaseModel):
    """Surface, delimited by curve loops."""

    curve_loop_tags: Sequence[StrictInt]
    tag: StrictInt

    class Config:
        allow_mutation = False


class PhysicalGroup(BaseModel):
    """Group of entities belonging to the same physical group."""

    dim: PhysicalDimension
    entity_tags: Sequence[StrictInt]
    tag: StrictInt

    class Config:
        allow_mutation = False


class SurfaceLoop(BaseModel):
    """Closed loop of surfaces."""

    surface_tags: Sequence[StrictInt]
    tag: StrictInt

    class Config:
        allow_mutation = False


class Volume(BaseModel):
    """Volume, delimited by surface loops."""

    surface_loop_tags: Sequence[StrictInt]
    tag: StrictInt

    class Config:
        allow_mutation = False


class GmshEntities(BaseModel):
    """Entities needed to create a mesh using gmsh."""

    points: Sequence[Point]
    lines: Sequence[Line]
    curve_loops: Sequence[CurveLoop]
    plane_surfaces: Optional[Sequence[PlaneSurface]]
    physical_groups: Optional[Sequence[PhysicalGroup]]


def mesh(
    gmsh_entities: GmshEntities,
    file_name: Path | str,
    dim: PhysicalDimension = 2,
    write_geo: bool = True,
    correct_curve_loops: bool = False,
    save_all: bool = False,
) -> None:
    """
    Create a mesh using gmsh.

    Args:
        gmsh_entities:
            Entities that represent the geometry.
        file_name:
            Path of the saved mesh.
        dim:
            Dimension of the mesh. Allowed values are 2 or 3.
        write_geo:
            Whether or not the .GEO file is saved too.
        correct_line_loops:
            Apparently gmsh offeres the possibility to correct
            line loops. TODO: Check if it is useful.
        save_all:
            If False only mesh elements are saved that belong to some
            physical group, as long there is at least one physical group.
    """
    assert dim in {2, 3}
    file_name = Path(file_name)
    gmsh.initialize()
    gmsh.model.add(file_name.name)

    for point in gmsh_entities.points:
        gmsh.model.geo.addPoint(
            x=point.x, y=point.y, z=point.z, meshSize=point.lc, tag=point.tag
        )

    for line in gmsh_entities.lines:
        gmsh.model.geo.addLine(
            startTag=line.start_tag, endTag=line.end_tag, tag=line.tag
        )

    for curve_loop in gmsh_entities.curve_loops:
        gmsh.model.geo.addCurveLoop(
            curveTags=curve_loop.line_tags,
            tag=curve_loop.tag,
            reorient=correct_curve_loops,
        )

    for plane_surface in gmsh_entities.plane_surfaces:
        gmsh.model.geo.addPlaneSurface(
            wireTags=plane_surface.curve_loop_tags, tag=plane_surface.tag
        )

    # Add physical groups
    for physical_group in gmsh_entities.physical_groups:
        gmsh.model.geo.addPhysicalGroup(
            dim=physical_group.dim.value,
            tag=physical_group.tag,
            tags=physical_group.entity_tags,
        )

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(dim)

    # Choose Dolfin compatible version
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    if save_all:
        # Write all mesh entities
        gmsh.option.setNumber("Mesh.SaveAll", 1)

    # Write the mesh file
    gmsh.write(str(file_name.with_suffix(".msh")))

    if write_geo:
        # Write the geo file
        gmsh.write(str(file_name.with_suffix(".geo_unrolled")))

    gmsh.finalize()


def write_geo(
    file_name: Path | str,
    gmsh_entities: GmshEntities,
    correct_curve_loops: bool = False,
):
    """Create a .GEO file."""
    gmsh.initialize()
    file_name = Path(file_name)
    gmsh.model.add(file_name.name)

    for point in gmsh_entities.points:
        gmsh.model.geo.addPoint(
            x=point.x, y=point.y, z=point.z, meshSize=point.lc, tag=point.tag
        )

    for line in gmsh_entities.lines:
        gmsh.model.geo.addLine(
            startTag=line.start_tag, endTag=line.end_tag, tag=line.tag
        )

    for curve_loop in gmsh_entities.curve_loops:
        gmsh.model.geo.addCurveLoop(
            curveTags=curve_loop.line_tags,
            tag=curve_loop.tag,
            reorient=correct_curve_loops,
        )

    for plane_surface in gmsh_entities.plane_surfaces:
        gmsh.model.geo.addPlaneSurface(
            wireTags=plane_surface.curve_loop_tags, tag=plane_surface.tag
        )

    for physical_group in gmsh_entities.physical_groups:
        gmsh.model.geo.addPhysicalGroup(
            dim=physical_group.dim.value,
            tag=physical_group.tag,
            tags=physical_group.entity_tags,
        )
    gmsh.model.geo.removeAllDuplicates()

    gmsh.model.geo.synchronize()

    # Write geo file
    gmsh.write(str(file_name.with_suffix(".geo_unrolled")))
    gmsh.finalize()
