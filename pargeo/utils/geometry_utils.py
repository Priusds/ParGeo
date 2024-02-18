import math
from typing import Callable

from shapely import box as shapely_box
from shapely import union_all
from shapely.affinity import translate as shapely_translate

from pargeo.utils.constants import DEFAULT_GRID_SIZE
from pargeo.utils.typing_utils import (
    GeometryCollection,
    MultiPolygon,
    Polygon,
    SubDomain,
    Vector,
)


def discretize_ellipse(
    midpoint: Vector, axis: Vector, angle: float, refs: int
) -> Polygon:
    """Discretize the ellipse.

    Args:
        midpoint: Ellipse's midpoint.
        axis: Ellipse's axis.
        angle: Rotation angle of the ellipse in radiant.
        refs: Number of discretization points.

    Return:
        List of boundary points in ccw order.
    """
    angles = [i * 2 * math.pi / refs for i in range(refs)]

    transformed_coords = rotate_counterclockwise(
        points=[polar_to_cartesian(angle, axis[0], axis[1]) for angle in angles],
        angle=angle,
    )
    mx, my = midpoint

    return Polygon([(x + mx, y + my) for x, y in transformed_coords])


def trigonometric_function(
    coefficients: list[Vector], scale: float
) -> Callable[[float], float]:
    """Returns the exponential of a trigonometric function.

    Args:
        coefficients: Coefficients [(a_i, b_i)_i] of the trigonometric function.
        scale: Scaling factor.

    The trigonometric function f is defined as:

        `f(x) = scale * exp(f_0(x) + ... + f_{n-1}(x))`

    where the basis functions f_i are defined as:

        `f_i(x) := a_i*cos(i*x) + b_i*sin(i*x)`
    """
    return lambda x: scale * math.exp(
        sum(
            a * math.cos(i * x) + b * math.sin(i * x)
            for i, (a, b) in enumerate(coefficients)
        )
    )


def rotate_counterclockwise(points: list[Vector], angle: float) -> list[Vector]:
    """Rotates points counter-clockwise around (0,0)."""
    cos_alpha = math.cos(angle)
    sin_alpha = math.sin(angle)
    return [
        (cos_alpha * x - sin_alpha * y, sin_alpha * x + cos_alpha * y)
        for x, y in points
    ]


def polar_to_cartesian(
    angle: float,
    rad_x: float,
    rad_y: float | None = None,
) -> Vector:
    """Coordinate transform, from polar to cartesian.

    Args:
        angles: List of angles.
        radii_x: List of radii, for the x-coordinate.
        radii_y: List of radii, for the y-coordinate.
    """
    if rad_y is None:
        rad_y = rad_x
    return (math.cos(angle) * rad_x, math.sin(angle) * rad_y)


def repeat(
    subdomain: SubDomain,
    repeat_info: tuple[Vector, float, int, bool],
    repeat_info_2: tuple[Vector, float, int, bool] | None = None,
    clipping_plane: SubDomain | None = None,
    grid_size: float = DEFAULT_GRID_SIZE,
) -> SubDomain:
    """Repeat the subdomain in one or two directions."""
    v_dir, v_scalar, v_rep, v_do_opposite = repeat_info
    if repeat_info_2 is not None:
        w_dir, w_scalar, w_rep, w_do_opposite = repeat_info_2
    else:
        w_dir, w_scalar, w_rep, w_do_opposite = (0, 0), 0, 0, False
    shifted_polygons = []

    bounds = subdomain.bounds
    i_range = range(v_rep) if not v_do_opposite else range(-v_rep + 1, v_rep)
    j_range = range(w_rep) if not w_do_opposite else range(-w_rep + 1, w_rep)
    for i in i_range:
        for j in j_range:
            x_off = i * v_scalar * v_dir[0]
            y_off = i * v_scalar * v_dir[1]
            x_off += j * w_scalar * w_dir[0]
            y_off += j * w_scalar * w_dir[1]

            # If clip to clipping_plane if given.
            if clipping_plane is not None:
                bounds_shifted = (
                    bounds[0] + x_off,
                    bounds[1] + y_off,
                    bounds[2] + x_off,
                    bounds[3] + y_off,
                )
                # Check if the shifted bounds intersect with the clipping plane.
                if not shapely_box(*bounds_shifted).intersects(clipping_plane):
                    continue

            shifted_polygons.append(
                shapely_translate(subdomain, xoff=x_off, yoff=y_off)
            )

    # Merge all shifted polygons.
    repeated_polygon = union_all(shifted_polygons)

    # Clip to domain.
    if clipping_plane is not None:
        repeated_polygon = repeated_polygon.intersection(
            clipping_plane, grid_size=grid_size
        )
        if isinstance(repeated_polygon, GeometryCollection):
            repeated_polygon = MultiPolygon(
                [poly for poly in repeated_polygon if poly.area > grid_size]
            )

    return repeated_polygon
