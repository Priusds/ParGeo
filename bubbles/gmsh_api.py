"""Module for mesh creation using gmsh."""
from enum import IntEnum
from pathlib import Path
from typing import Optional, Sequence
from itertools import groupby

import gmsh
from pydantic import BaseModel, StrictInt
import shapely

from bubbles.two_d.topology import Bubble
from bubbles.two_d.topology_new import Topology


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


def bubble_to_gmsh_entities(
    polygon, level, point_tag, line_tag, curve_loop_tag, plane_surface_tag
):
    #assert not bubble.is_hole
    points_total = []
    lines_total = []
    curve_loops = []

    physical_group_tag = level + 1  # +1 because physical group tags start at 1
    # ==================================================
    # Write exterior loop
    # ==================================================
    assert isinstance(polygon, shapely.Polygon)
    # Add points
    points_local = []
    for x, y in polygon.exterior.coords[:-1]:
        point = Point(x=x, y=y, z=0.0, lc=0.1, tag=point_tag)
        points_local.append(point)
        point_tag += 1

    # Add lines
    lines_local = []
    for i in range(len(points_local) - 1):
        line = Line(
            start_tag=points_local[i].tag, end_tag=points_local[i + 1].tag, tag=line_tag
        )
        lines_local.append(line)
        line_tag += 1
    line = Line(
        start_tag=points_local[-1].tag, end_tag=points_local[0].tag, tag=line_tag
    )
    lines_local.append(line)
    line_tag += 1

    points_total = points_total + points_local
    lines_total = lines_total + lines_local
    # Add curve loop
    curve_loop = CurveLoop(
        line_tags=[line.tag for line in lines_local], tag=curve_loop_tag
    )
    curve_loops.append(curve_loop)
    curve_loop_tag += 1

    # ==================================================
    # Write interior loops
    # ==================================================
    for loop in polygon.interiors:
        points_local = []
        lines_local = []
        for x, y in loop.coords[:-1]:
            point = Point(x=x, y=y, z=0.0, lc=0.1, tag=point_tag)
            points_local.append(point)
            point_tag += 1

        # Add lines
        for i in range(len(points_local) - 1):
            line = Line(
                start_tag=points_local[i].tag,
                end_tag=points_local[i + 1].tag,
                tag=line_tag,
            )
            lines_local.append(line)
            line_tag += 1
        line = Line(
            start_tag=points_local[-1].tag, end_tag=points_local[0].tag, tag=line_tag
        )
        lines_local.append(line)
        line_tag += 1

        # Add curve loop
        curve_loop = CurveLoop(
            line_tags=[line.tag for line in lines_local], tag=curve_loop_tag
        )
        curve_loops.append(curve_loop)
        curve_loop_tag += 1

        points_total = points_total + points_local
        lines_total = lines_total + lines_local

    # ==================================================
    # Write surface
    # ==================================================
    plane_surfaces = [
        PlaneSurface(
            curve_loop_tags=[loop.tag for loop in curve_loops], tag=plane_surface_tag
        )
    ]
    physical_groups = [
        PhysicalGroup(dim=2, entity_tags=[plane_surface_tag], tag=physical_group_tag)
    ]
    plane_surface_tag += 1

    return (
        GmshEntities(
            points=points_total,
            lines=lines_total,
            curve_loops=curve_loops,
            plane_surfaces=plane_surfaces,
            physical_groups=physical_groups,
        ),
        point_tag,
        line_tag,
        curve_loop_tag,
        plane_surface_tag,
    )


def topology_to_gmsh_entities(topo: Topology):
    all_entities_list = []

    bubbles = topo.flatten()

    point_tag = 1
    line_tag = 1
    curve_loop_tag = 1
    plane_surface_tag = 1

    for polygon, level in bubbles:
        if not level in topo.hole_levels:
            (
                gmsh_entities,
                point_tag,
                line_tag,
                curve_loop_tag,
                plane_surface_tag,
            ) = bubble_to_gmsh_entities(
                polygon, level, point_tag, line_tag, curve_loop_tag, plane_surface_tag
            )
            all_entities_list.append(gmsh_entities)

    # Group all physical groups by tag
    physical_groups_merged = []
    physical_groups = [e.physical_groups[0] for e in all_entities_list]
    physical_groups.sort(key=lambda x: x.tag)

    for tag, physical_group in groupby(physical_groups, lambda x: x.tag):
        entity_tags = [s.entity_tags[0] for s in physical_group]
        physical_groups_merged.append(
            PhysicalGroup(dim=2, entity_tags=entity_tags, tag=tag)
        )

    return GmshEntities(
        points=[i for ent in all_entities_list for i in ent.points],
        lines=[i for ent in all_entities_list for i in ent.lines],
        curve_loops=[i for ent in all_entities_list for i in ent.curve_loops],
        plane_surfaces=[i for ent in all_entities_list for i in ent.plane_surfaces],
        physical_groups=physical_groups_merged,
    )


def bubbles_to_gmsh_entities(bubbles: list[Bubble]):
    all_entities_list = []

    point_tag = 1
    line_tag = 1
    curve_loop_tag = 1
    plane_surface_tag = 1

    for bubble in bubbles:
        if not bubble.is_hole:
            (
                gmsh_entities,
                point_tag,
                line_tag,
                curve_loop_tag,
                plane_surface_tag,
            ) = bubble_to_gmsh_entities(
                polygon, level, point_tag, line_tag, curve_loop_tag, plane_surface_tag
            )
            all_entities_list.append(gmsh_entities)

    # Group all physical groups by tag
    physical_groups_merged = []
    physical_groups = [e.physical_groups[0] for e in all_entities_list]
    physical_groups.sort(key=lambda x: x.tag)

    for tag, physical_group in groupby(physical_groups, lambda x: x.tag):
        entity_tags = [s.entity_tags[0] for s in physical_group]
        physical_groups_merged.append(
            PhysicalGroup(dim=2, entity_tags=entity_tags, tag=tag)
        )

    return GmshEntities(
        points=[i for ent in all_entities_list for i in ent.points],
        lines=[i for ent in all_entities_list for i in ent.lines],
        curve_loops=[i for ent in all_entities_list for i in ent.curve_loops],
        plane_surfaces=[i for ent in all_entities_list for i in ent.plane_surfaces],
        physical_groups=physical_groups_merged,
    )
