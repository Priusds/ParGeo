"""API for Gmsh.

This module provides an API for Gmsh. It allows the user to write a Gmsh file 
and create a mesh using Gmsh.
"""

from itertools import groupby
from pathlib import Path
from typing import NamedTuple, Sequence, TypedDict

import gmsh
from shapely import Polygon

from pargeo.domain import Domain

PHISICAL_DIMENSION = 2


class Point(NamedTuple):
    """Two dimensional point."""

    x: float
    y: float
    z: float
    lc: float
    tag: int


class Line(NamedTuple):
    """Line between two points."""

    start_tag: int
    end_tag: int
    tag: int


class CurveLoop(NamedTuple):
    """Closed loop of lines."""

    line_tags: Sequence[int]
    tag: int


class PlaneSurface(NamedTuple):
    """Surface, delimited by curve loops."""

    curve_loop_tags: Sequence[int]
    tag: int


class PhysicalGroup(NamedTuple):
    """Group of entities belonging to the same physical group."""

    entity_tags: Sequence[int]
    tag: int


class SurfaceLoop(NamedTuple):
    """Closed loop of surfaces."""

    surface_tags: Sequence[int]
    tag: int


class Volume(NamedTuple):
    """Volume, delimited by surface loops."""

    surface_loop_tags: Sequence[int]
    tag: int


class GmshEntities(TypedDict):
    """Entities needed to create a mesh using gmsh."""

    points: Sequence[Point]
    lines: Sequence[Line]
    curve_loops: Sequence[CurveLoop]
    plane_surfaces: Sequence[PlaneSurface]
    physical_groups: Sequence[PhysicalGroup]
    surface_loops: Sequence[SurfaceLoop]
    volumes: Sequence[Volume]


def write_msh_from_entities(
    gmsh_entities: GmshEntities,
    file_name: Path,
    write_geo: bool = True,
    correct_curve_loops: bool = False,
    save_all: bool = False,
) -> None:
    """
    Write a mesh .MSH file using GMSH.

    Args:
        gmsh_entities: Entities that represent the geometry.
        file_name: Path of the saved mesh, the extension is added automatically.
        write_geo: Whether or not the .GEO file is saved too.
        correct_line_loops: Apparently gmsh offeres the possibility to correct
            line loops. TODO: Check if it is useful.
        save_all: If False only mesh elements are saved that belong to some
            physical group, as long there is at least one physical group.
    """

    gmsh.initialize()
    gmsh.model.add(file_name.name)

    for point in gmsh_entities["points"]:
        gmsh.model.geo.addPoint(
            x=point.x, y=point.y, z=point.z, meshSize=point.lc, tag=point.tag
        )

    for line in gmsh_entities["lines"]:
        gmsh.model.geo.addLine(
            startTag=line.start_tag, endTag=line.end_tag, tag=line.tag
        )

    for curve_loop in gmsh_entities["curve_loops"]:
        gmsh.model.geo.addCurveLoop(
            curveTags=curve_loop.line_tags,
            tag=curve_loop.tag,
            reorient=correct_curve_loops,
        )

    for plane_surface in gmsh_entities["plane_surfaces"]:
        gmsh.model.geo.addPlaneSurface(
            wireTags=plane_surface.curve_loop_tags, tag=plane_surface.tag
        )

    # Add physical groups
    for physical_group in gmsh_entities["physical_groups"]:
        gmsh.model.geo.addPhysicalGroup(
            dim=PHISICAL_DIMENSION,
            tag=physical_group.tag,
            tags=physical_group.entity_tags,
        )

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(PHISICAL_DIMENSION)

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


def write_geo_from_entities(
    file_name: Path,
    gmsh_entities: GmshEntities,
    correct_curve_loops: bool = False,
):
    """Write a file readable by GMSH.

    The geo file can be used in gmsh.

    Args:
        file_name: Name of the file. The extension ".geo_unrolled" is added
            automatically.
        gmsh_entities: Entities that represent the geometry.
        correct_curve_loops: Apparently gmsh offeres the possibility to correct
            line loops. TODO: Check if it is useful.
    """
    gmsh.initialize()
    gmsh.model.add(file_name.name)

    for point in gmsh_entities["points"]:
        gmsh.model.geo.addPoint(
            x=point.x, y=point.y, z=point.z, meshSize=point.lc, tag=point.tag
        )

    for line in gmsh_entities["lines"]:
        gmsh.model.geo.addLine(
            startTag=line.start_tag, endTag=line.end_tag, tag=line.tag
        )

    for curve_loop in gmsh_entities["curve_loops"]:
        gmsh.model.geo.addCurveLoop(
            curveTags=curve_loop.line_tags,
            tag=curve_loop.tag,
            reorient=correct_curve_loops,
        )

    for plane_surface in gmsh_entities["plane_surfaces"]:
        gmsh.model.geo.addPlaneSurface(
            wireTags=plane_surface.curve_loop_tags, tag=plane_surface.tag
        )

    for physical_group in gmsh_entities["physical_groups"]:
        gmsh.model.geo.addPhysicalGroup(
            dim=PHISICAL_DIMENSION,
            tag=physical_group.tag,
            tags=physical_group.entity_tags,
        )
    gmsh.model.geo.removeAllDuplicates()

    gmsh.model.geo.synchronize()

    # Write geo file
    gmsh.write(str(file_name.with_suffix(".geo_unrolled")))
    gmsh.finalize()


def polygon_to_entities(
    polygon: Polygon,
    level: int,
    point_tag: int,
    line_tag: int,
    curve_loop_tag: int,
    plane_surface_tag: int,
):
    """Convert a shapely polygon to gmsh entities.

    Args:
        polygon: Shapely polygon.
        level: Physical group of the surface.
        point_tag: Lower bound for free point tags.
        line_tag: Lower bound for free line tags.
        curve_loop_tag: Lower bound for free curve_loop tags.
        plane_surface_tag: Lower bound for free surface tags.
    """
    # assert not bubble.is_hole
    points_total: list[Point] = []
    lines_total: list[Line] = []
    curve_loops: list[CurveLoop] = []

    physical_group_tag = level + 1  # +1 because physical group tags start at 1
    # ==================================================
    # Write exterior loop
    # ==================================================
    assert isinstance(polygon, Polygon)
    # Add points
    points_local: list[Point] = []
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

    points_total.extend(points_local)
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
        PhysicalGroup(
            entity_tags=[plane_surface_tag],
            tag=physical_group_tag,
        )
    ]
    plane_surface_tag += 1

    entities: GmshEntities = {
        "points": points_total,
        "lines": lines_total,
        "curve_loops": curve_loops,
        "plane_surfaces": plane_surfaces,
        "physical_groups": physical_groups,
        "surface_loops": [],
        "volumes": [],
    }
    return (
        entities,
        point_tag,
        line_tag,
        curve_loop_tag,
        plane_surface_tag,
    )


def domain_to_entities(domain: Domain) -> GmshEntities:
    """Convert a domain to gmsh entities."""
    all_entities_list = []

    # TODO: Rename this to something more meaningful
    pargeo = domain.as_list()

    point_tag = 1
    line_tag = 1
    curve_loop_tag = 1
    plane_surface_tag = 1

    for polygon, level in pargeo:
        if level not in domain.holes:
            (
                gmsh_entities,
                point_tag,
                line_tag,
                curve_loop_tag,
                plane_surface_tag,
            ) = polygon_to_entities(
                polygon, level, point_tag, line_tag, curve_loop_tag, plane_surface_tag
            )
            all_entities_list.append(gmsh_entities)

    # Group all physical groups by tag
    physical_groups_merged = []
    physical_groups = [e["physical_groups"][0] for e in all_entities_list]
    physical_groups.sort(key=lambda x: x.tag)

    for tag, physical_group in groupby(physical_groups, lambda x: x.tag):
        entity_tags = [s.entity_tags[0] for s in physical_group]
        physical_groups_merged.append(PhysicalGroup(entity_tags=entity_tags, tag=tag))

    gmsh_entities_all: GmshEntities = {
        "points": [i for ent in all_entities_list for i in ent["points"]],
        "lines": [i for ent in all_entities_list for i in ent["lines"]],
        "curve_loops": [i for ent in all_entities_list for i in ent["curve_loops"]],
        "plane_surfaces": [
            i for ent in all_entities_list for i in ent["plane_surfaces"]
        ],
        "physical_groups": physical_groups_merged,
        "surface_loops": [],
        "volumes": [],
    }
    return gmsh_entities_all


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
