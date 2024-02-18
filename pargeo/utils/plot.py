from typing import Mapping

import matplotlib.pyplot as plt

from pargeo.utils.constants import BACKGROUND_LEVEL
from pargeo.utils.typing_utils import Color, Level, MultiPolygon, Polygon, SubDomain


def plot_legend(
    holes: set[int],
    colormap: Mapping[Level, Color],
    background_level: Level = BACKGROUND_LEVEL,
):
    """Make a legend."""
    handles = []
    descriptions = []
    for lvl, color in colormap.items():
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, edgecolor="black"))
        if lvl in holes:
            descriptions.append(f"{lvl}. level (hole)")
        else:
            if lvl == background_level:
                descriptions.append(f"Background (level = {background_level})")
            else:
                descriptions.append(f"{lvl}. level")
    plt.legend(
        handles,
        descriptions,
        loc="upper left",
        bbox_to_anchor=(0.95, 1),
    ).set_draggable(True)


def plot_polygon(
    polygon: Polygon,
    color: Color,
    boundary_color: Color,
    boundary_style: str,
    linewidth: float,
    show: bool = False,
):
    """Plot a polygon."""
    x, y = polygon.exterior.xy
    plt.fill(x, y, color=color)

    # plot the exterior boundary
    if not boundary_color == "none":
        plt.plot(x, y, boundary_style, linewidth=linewidth, color=boundary_color)

    # plot interiors in white
    for interior in polygon.interiors:
        x, y = interior.xy
        plt.fill(x, y, color="white")
        plt.plot(x, y, "-", linewidth=1, color="black")
    if show:
        plt.show()


def plot_subdomain(subdomain: SubDomain):
    color = "blue"
    boundary_color = "black"
    boundary_style = "-"
    linewidth = 1.0
    show = True
    if isinstance(subdomain, MultiPolygon):
        for polygon in subdomain.geoms:
            plot_polygon(polygon, color, boundary_color, boundary_style, linewidth)
    else:
        plot_polygon(subdomain, color, boundary_color, boundary_style, linewidth)
    if show:
        plt.show()
