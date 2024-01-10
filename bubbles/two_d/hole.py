"""Define two-dimensional bubbles here.

TODO: Refactor this module. Do following:
    - Rename `Hole`
    - Integrate Shapely
    - Rename `Rectangular` to `Rectangle` and add more boundary points.
    - Maybe move Bubble class here.
"""
from abc import ABC, abstractmethod
from typing import Any, Callable, Optional

import numpy as np
from numpy.typing import ArrayLike


class Hole(ABC):
    """Abstract class for two-dimensional bubbles."""

    def __init__(self, type_):
        self.type = type_

    @abstractmethod
    def discretize_hole(self, refs: int) -> list[tuple[float, float]]:
        """Discretize the boundary of the hole.

        Return the boundary points in counter-clockwise order.

        Args:
            refs:
                Number of boundary points that are returned.
        """
        raise NotImplementedError

    @abstractmethod
    def to_dict(self) -> dict:
        """Returns the hole information in dict format."""
        raise NotImplementedError

    def info(self):
        """Prints information about the hole."""
        print(self.to_dict())

    def get_point(self, angle: float) -> tuple[float, float]:
        """Returns the corresponding boundary point.

        Note: So far all holes are stellar-like polygons, which means
        that every hole has a midpoint, and to a given angle corresponds
        a boundary point, which is returned by this method.
        """
        raise NotImplementedError


class Circle(Hole):
    """Circle bubble."""

    def __init__(self, midpoint: tuple[float, float], radius: float):
        """Creates a circle.

        Args:
            midpoint:
                The midpoint.
            radius:
                Length of the radius.
        """
        super().__init__("circle")
        self.midpoint = midpoint
        self.radius = radius

    def to_dict(self):
        """Returns the hole information in dict format."""
        return {"type": self.type, "midpoint": self.midpoint, "radius": self.radius}

    def discretize_hole(self, refs: int) -> list[tuple[float, float]]:
        """Discretize the boundary of the hole.

        Return the boundary points in counter-clockwise order.

        Args:
            refs:
                Number of boundary points that are returned.
        """
        return discretize_ellipse(self.midpoint, (self.radius, self.radius), 0, refs)

    def get_point(self, angle, axis):
        """Returns the corresponding boundary point."""
        if axis == 0:
            return self.radius * np.cos(angle) + self.midpoint[0]
        elif axis == 1:
            return self.radius * np.sin(angle) + self.midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Ellipse(Hole):
    """Ellipse bubble."""

    def __init__(
        self, midpoint: tuple[float, float], axis: tuple[float, float], angle: float
    ):
        """Creates an ellipse.

        Args:
            midpoint:
                The midpoint.
            axis:
                Length of the two axis.
            angle:
                In radian, rotates the ellipse counter-clockwise.
        """
        super().__init__("ellipse")
        self.midpoint = midpoint
        self.axis = axis
        self.angle = angle

    def to_dict(self):
        """Returns the hole information in dict format."""
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "axis": self.axis,
            "angle": self.angle,
        }

    def discretize_hole(self, refs: int) -> list[tuple[float, float]]:
        """Discretize the boundary of the hole.

        Return the boundary points in counter-clockwise order.

        Args:
            refs:
                Number of boundary points that are returned.
        """
        return discretize_ellipse(self.midpoint, self.axis, self.angle, refs)

    def get_point(self, angle, axis):
        """Returns the corresponding boundary point."""
        if axis == 0:
            return self.axis[0] * np.cos(angle) + self.midpoint[0]
        elif axis == 1:
            return self.axis[1] * np.sin(angle) + self.midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Stellar(Hole):
    """Stellar bubble."""

    def __init__(
        self,
        midpoint: tuple[float, float],
        radius: float,
        coefficient: Optional[np.ndarray] = None,
    ):
        """Creates a star-like geometry.

        Args:
            midpoint:
                The midpoint.
            radius:
                Approximately the average radius.
            coefficient: numpy array of shape (n, 2) with n > 0.
                Hint: for smooth stellar holes, you can çhoose k from {.1, .2}
                and then `coeffients = np.random.uniform(-k, k, (int(1/k), 2))`
        """
        super().__init__("stellar")
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

    def to_dict(self):
        """Returns the hole information in dict format."""
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "coefficient": self.coefficient.tolist(),
            "radius": self.radius,
        }

    def discretize_hole(self, refs: int) -> list[tuple[float, float]]:
        """Discretize the boundary of the hole.

        Return the boundary points in counter-clockwise order.

        Args:
            refs:
                Number of boundary points that are returned.
        """
        return discetize_stellar_polygon(
            self.midpoint, self.radius, self.coefficient, refs
        )

    def get_point(self, angle, axis):
        """Returns the corresponding boundary point."""
        r = self.f(angle)

        if axis == 0:
            return r * np.cos(angle) + self.midpoint[0]
        elif axis == 1:
            return r * np.sin(angle) + self.midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Rectangular(Hole):
    """Rectangular bubble."""

    def __init__(self, midpoint: tuple[float, float], width: float, height: float):
        """Creates a rectangle.

        Args:
            midpoint:
                The midpoint.
            width:
                Width of the rectangle.
            height:
                Height of the rectangle.
        """
        super().__init__("rectangular")
        self.midpoint = midpoint
        self.width = width
        self.height = height

    def to_dict(self):
        """Returns the hole information in dict format."""
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "width": self.width,
            "height": self.height,
        }

    def discretize_hole(self, refs: Any) -> list[tuple[float, float]]:
        """Discretize the boundary of the hole.

        Return only the corner points in counter-clockwise order.

        Args:
            refs:
                Not used.
        """
        midpoint = self.midpoint
        width = self.width
        height = self.height
        p0 = (midpoint[0] - width / 2, midpoint[1] - height / 2)
        p1 = (midpoint[0] + width / 2, midpoint[1] - height / 2)
        p2 = (midpoint[0] + width / 2, midpoint[1] + height / 2)
        p3 = (midpoint[0] - width / 2, midpoint[1] + height / 2)

        return [p0, p1, p2, p3]


def discretize_ellipse(
    midpoint: tuple[float, float], axis: tuple[float, float], angle: float, refs: int
) -> list[tuple[float, float]]:
    """Discretize the ellipse.

    Args:
        midpoint:
            Ellipse's midpoint.
        axis:
            Ellipse's axis.
        angle:
            Rotation angle of the ellipse in radiant.
        refs:
            Number of discretization points.

    Return:
        List of boundary points in counter-clockwise order.
    """
    angles = np.linspace(0, 2 * np.pi, refs, endpoint=False)
    coords = polar_to_cartesian(angles, axis[0], axis[1])
    transformed_coords = rotate_counterclockwise(coords, angle) + midpoint

    return [tuple(coords) for coords in transformed_coords]  # type: ignore


def discretize_ellipse_noisy(
    midpoint: tuple[float, float],
    axis: tuple[float, float],
    angle: float,
    refs: int,
    var: float,
) -> list[tuple[float, float]]:
    """Discretize the ellipse with additional noise.

    Args:
        midpoint:
            Ellipse's midpoint.
        axis:
            Ellipse's axis.
        angle:
            Rotation angle of the ellipse in radiant.
        refs:
            Number of discretization points.
        var:
            Standard deviation of the Gaussian noise.
            Hint: choose var ~ axis/10
    Return:
        List of boundary points in counter-clockwise order.
    """
    noise = np.random.normal(0, var, refs)

    angles = np.linspace(0, 2 * np.pi, refs, endpoint=False)
    coords = polar_to_cartesian(angles, axis[0] + noise, axis[1] + noise)
    transformed_coords = rotate_counterclockwise(coords, angle) + midpoint

    return [tuple(coords) for coords in transformed_coords]  # type: ignore


def discetize_stellar_polygon(
    midpoint: tuple[float, float], radius: float, coefficients: np.ndarray, refs: int
) -> list[tuple[float, float]]:
    """Discretize the stellar polygon.

    Args:
        midpoint:
            The midpoint.
        radius:
            The radius.
        coefficients:
            The coefficients.
        refs:
            Number of discretization points.
    """
    angles = np.linspace(0, 2 * np.pi, refs, endpoint=False)
    f = trigonometric_function(coefficients, radius)
    polygon = polar_to_cartesian(angles, [f(angle) for angle in angles]) + np.array(
        midpoint
    )

    return [tuple(coords) for coords in polygon]  # type: ignore


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
        scale:
            Scaling factor.
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


def rotate_counterclockwise(p, alpha):
    """Rotates points counter-clockwise around the origin."""
    p = np.array(p)
    q = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
    return p @ q.transpose()


def rotate_clockwise(p, alpha):
    """Rotates points clockwise around the origin."""
    # TODO: Check if needed.
    p = np.array(p)
    q_inv = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    return q_inv @ p


def polar_to_cartesian(
    angles: ArrayLike,
    radii_x: ArrayLike,
    radii_y: Optional[ArrayLike] = None,
) -> np.ndarray:
    """Coordinate transform, from polar to cartesian.

    x_i, y_i = cos(angle_i)*rx_i, sin(angle_i)*ry_i

    Args:
        angles: array-like
            List of angles.
        radii_x: array-like
            List of radii, for the x-coordinate.
        radii_y: array-like
            List of radii, for the y-coordinate.
    """
    if radii_y is None:
        radii_y = radii_x
    return np.vstack((np.cos(angles) * radii_x, np.sin(angles) * radii_y)).transpose()
