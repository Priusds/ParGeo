from typing import Mapping, Sequence, Union

import matplotlib.pyplot as plt

from pargeo.utils.constants import BACKGROUND_LEVEL
from pargeo.utils.typing_utils import Level, MultiPolygon, Polygon, SubDomain
from matplotlib.colors import Normalize

Color = Union[str, tuple[int, int, int], tuple[int, int, int, int]]


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


def plot_subdomain(
    subdomain: SubDomain,
    color="blue",
    boundary_color="black",
    boundary_style="-",
    linewidth=1.0,
    show: bool = True,
):
    """Plot a subdomain."""
    if isinstance(subdomain, MultiPolygon):
        for polygon in subdomain.geoms:
            plot_polygon(polygon, color, boundary_color, boundary_style, linewidth)
    else:
        plot_polygon(subdomain, color, boundary_color, boundary_style, linewidth)
    if show:
        plt.show()


class DefaultColors:
    """Default colors for plotting."""

    boundary = "black"
    domain_boundary = "black"
    domain = "silver"
    hole = "white"
    interior_filling = "white"
    inter_hole_bound = "white"

    @staticmethod
    def get_color_map(
        levels: Sequence[int],
        holes: set[int],
        color_holes=False,
        background_level: int = BACKGROUND_LEVEL,
    ) -> Mapping[Level, Color]:
        """Get a color map for the levels."""
        lvl2cl = dict()
        level_set = set(levels)
        if background_level in level_set:
            lvl2cl[background_level] = DefaultColors.domain
            level_set.remove(background_level)
        if len(level_set) > 0:
            norm = Normalize(vmin=min(level_set), vmax=max(level_set))
            for lvl in level_set:
                if lvl in holes and not color_holes:
                    lvl2cl[lvl] = DefaultColors.hole
                else:
                    lvl2cl[lvl] = plt.cm.cool(norm(lvl))  # type: ignore
        return lvl2cl
