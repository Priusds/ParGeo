"""Module containing utility methods."""
from typing import Any

from bubbles.gmsh_api import (
    CurveLoop,
    Line,
    PlaneSurface,
    Point,
    SurfaceLoop,
    Volume,
)


def geometry_to_gmsh_entities(geometry, default_lc=1) -> dict[str, Any]:
    """Transform old geometry structure in gmsh entity structure."""
    entities_dict = {
        "points": [],
        "lines": [],
        "curve_loops": [],
        "plane_surfaces": [],
        "surface_loops": [],
        "volumes": [],
        "physical_groups": [],
    }
    entities_dict["points"] = [
        Point(
            x=point["coords"][0],
            y=point["coords"][1],
            z=point["coords"][2],
            lc=point["lc"] if point["lc"] else default_lc,
            tag=id,
        )
        for id, point in geometry["points"].items()
    ]

    if "lines" in geometry.keys():
        entities_dict["lines"] = [
            Line(start_tag=start, end_tag=end, tag=id)
            for id, (start, end) in geometry["lines"].items()
        ]

    if "lineLoops" in geometry.keys():
        entities_dict["curve_loops"] = [
            CurveLoop(line_tags=line_ids, tag=id)
            for id, line_ids in geometry["lineLoops"].items()
        ]

    if "surfaces" in geometry.keys():
        entities_dict["plane_surfaces"] = [
            PlaneSurface(curve_loop_tags=curve_loop_ids, tag=id)
            for id, curve_loop_ids in geometry["surfaces"].items()
        ]

    if "surfaceLoops" in geometry.keys():
        entities_dict["surface_loops"] = [
            SurfaceLoop(surface_tags=surface_tags, tag=id)
            for id, surface_tags in geometry["surfaceLoops"].items()
        ]

    if "volumes" in geometry.keys():
        entities_dict["volumes"] = [
            Volume(surface_loop_tags=surface_loop_tags, tag=id)
            for id, surface_loop_tags in geometry["volumes"].items()
        ]

    if "physical_groups" in geometry.keys():
        raise NotImplementedError("Physical groupd not supported yet.")

    return entities_dict
