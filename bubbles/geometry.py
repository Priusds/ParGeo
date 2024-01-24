"""Geometry module.

TODO: Maybe Rename to `continuous_geometries.py`?

Create geometries that can be added to the topology.
"""
import math
from abc import ABC, abstractmethod
from typing import Callable, Optional

import numpy as np
from numpy.typing import ArrayLike
from shapely.geometry import Polygon


class Geometry(ABC):
    """Abstract class for geometries."""

    @abstractmethod
    def to_polygon(self, *Args) -> Polygon:
        """Return a shapely polygon representation of self."""
        pass


class Rectangle(Geometry):
    """Rectangular geometry."""

    def __init__(self, midpoint: tuple[float, float], width: float, height: float):
        """Creates a rectangle.

        Args:
            midpoint: The midpoint.
            width: Width of the rectangle.
            height: Height of the rectangle.
        """
        self.midpoint = midpoint
        self.width = width
        self.height = height

    def to_polygon(self) -> Polygon:
        """Return the corners of the rectangle in counter-clockwise order."""
        midpoint = self.midpoint
        width = self.width
        height = self.height
        p0 = (midpoint[0] - width / 2, midpoint[1] - height / 2)
        p1 = (midpoint[0] + width / 2, midpoint[1] - height / 2)
        p2 = (midpoint[0] + width / 2, midpoint[1] + height / 2)
        p3 = (midpoint[0] - width / 2, midpoint[1] + height / 2)

        return Polygon([p0, p1, p2, p3])


class StarLike(Geometry):
    def __init__(self, midpoint: tuple[float, float]) -> None:
        self.midpoint = midpoint

    @abstractmethod
    def distance(self, angle: float) -> float:
        """Returns the distance at a given angle to the midpoint."""
        raise NotImplementedError

    def to_polygon(self, refs) -> Polygon:
        """Return a shapely polygon representation of self."""
        return self.discretize(refs)

    def discretize(self, refs) -> Polygon:
        angles = np.linspace(0, 2 * np.pi, refs, endpoint=False)
        radii = [self.distance(angle) for angle in angles]
        x_coords = np.cos(angles) * radii + self.midpoint[0]
        y_coords = np.sin(angles) * radii + self.midpoint[1]

        return Polygon(zip(x_coords, y_coords))


class Circle(StarLike):
    def __init__(self, midpoint: tuple[float, float], radius) -> None:
        self.radius = radius
        super().__init__(midpoint)

    def distance(self, angle: float) -> float:
        return self.radius


class Ellipse(StarLike):
    def __init__(
        self, midpoint: tuple[float, float], axis: tuple[float, float], angle: float
    ) -> None:
        self.axis = axis
        self.angle = angle
        super().__init__(midpoint)

    def distance(self, angle: float) -> float:
        # TODO: Check if this is correct
        return math.sqrt(
            (self.axis[0] * math.cos(angle)) ** 2
            + (self.axis[1] * math.sin(angle)) ** 2
        )

    def discretize(self, refs) -> Polygon:
        return discretize_ellipse(self.midpoint, self.axis, self.angle, refs)


class Stellar(StarLike):
    """Stellar bubble."""

    def __init__(
        self,
        midpoint: tuple[float, float],
        radius: float,
        coefficient: Optional[np.ndarray] = None,
    ):
        """Creates a star-like geometry.

        Args:
            midpoint: The midpoint.
            radius: Approximately the average radius.
            coefficient: numpy array of shape (n, 2) with n > 0.
                Hint: for smooth stellar holes, you can çhoose k from {.1, .2}
                and then `coeffients = np.random.uniform(-k, k, (int(1/k), 2))`
        """
        self.midpoint = midpoint
        self.radius = radius
        if coefficient is None:
            k = 0.2
            j = int(1 / k)
            coefficient = np.random.uniform(-k, k, (j, 2))
        if len(coefficient) == 0:
            raise ValueError("`coefficient` length must be greater then zero.")
        self.coefficient = coefficient
        self.f = trigonometric_function(self.coefficient, self.radius)

    def distance(self, angle: float) -> float:
        return self.f(angle)


def discretize_ellipse(
    midpoint: tuple[float, float], axis: tuple[float, float], angle: float, refs: int
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
    angles = np.linspace(0, 2 * np.pi, refs, endpoint=False)
    coords = polar_to_cartesian(angles, axis[0], axis[1])
    transformed_coords = rotate_counterclockwise(coords, angle) + midpoint

    return Polygon([tuple(coords) for coords in transformed_coords])


def trigonometric_function(
    coefficients: np.ndarray, scale: float
) -> Callable[[float], float]:
    """Returns the exponential of a trigonometric function.

    The trigonometric function f is defined as:

        `f(x) = scale * exp(mean(a_0, b_0) + f_1(x) + ... + f_{n-1}(x))`

    where the basis functions f_i are defined as:

        `f_i(x) := a_i*cos(i*x) + b_i*sin(i*x)`

    Args:
        coefficients: numpy array of shape (n,2)
            Each row represents a pair (a_i, b_i).
        scale: Scaling factor.
    """
    if len(coefficients) == 0:
        raise ValueError("Length of `coefficients` must be greater then 0.")

    if len(coefficients) == 1:
        return lambda x: scale * np.exp(coefficients[0].mean())  # type: ignore

    else:
        i_matrix = np.array(range(1, len(coefficients)))
        return lambda x: scale * np.exp(  # type: ignore
            coefficients[0].mean()
            + np.dot(coefficients[1:, 0], np.cos(i_matrix * x))
            + np.dot(coefficients[1:, 1], np.sin(i_matrix * x))
        )


def rotate_counterclockwise(p: np.ndarray, alpha: float) -> np.ndarray:
    """Rotates points counter-clockwise around (0,0)."""
    p = np.array(p)
    q = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    return p @ q


def polar_to_cartesian(
    angles: ArrayLike,
    radii_x: ArrayLike,
    radii_y: Optional[ArrayLike] = None,
) -> np.ndarray:
    """Coordinate transform, from polar to cartesian.

    x_i, y_i = cos(angle_i)*rx_i, sin(angle_i)*ry_i

    Args:
        angles: List of angles.
        radii_x: List of radii, for the x-coordinate.
        radii_y: List of radii, for the y-coordinate.
    """
    if radii_y is None:
        radii_y = radii_x
    return np.vstack((np.cos(angles) * radii_x, np.sin(angles) * radii_y)).transpose()
