"""Geometry module.

Create geometries that can be added to the topology.
"""
import math
import random
from abc import ABC, abstractmethod
from typing import Callable

from shapely.geometry import Polygon


class Geometry(ABC):
    """Abstract class for geometries."""

    @abstractmethod
    def to_polygon(self, *args, **kwargs) -> Polygon:
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


class NStar(Geometry):
    """NStar geometry."""

    def __init__(
        self,
        midpoint: tuple[float, float],
        radius_in: float,
        radius_out: float,
        n_peaks: int,
        alpha: float = 0.0,
    ):
        """Make an NStar.

        Args:
            midpoint: The midpoint.
            radius_in: The inner radius.
            radius_out: The outer radius.
            N: Number of points.
            alpha: Angle offset. TODO: Explain this.
        """
        if radius_in > radius_out:
            raise ValueError(
                "The parameter setting radius_in > radius_out not allowed."
            )
        self._midpoint = midpoint
        self._inner_radius = radius_in
        self._outer_radius = radius_out
        self._n_peaks = n_peaks
        self.__alpha = alpha

    def to_polygon(self) -> Polygon:
        aux_angles = [2 * math.pi * i / self._n_peaks for i in range(self._n_peaks)]
        offset = math.pi / self._n_peaks

        angles_in = [angle + self.__alpha + offset for angle in aux_angles]
        angles_out = [angle + self.__alpha for angle in aux_angles]

        coords_in = [
            polar_to_cartesian(angle, self._inner_radius) for angle in angles_in
        ]
        coords_out = [
            polar_to_cartesian(angle, self._outer_radius) for angle in angles_out
        ]

        # Translate to midpoint
        mx, my = self._midpoint
        coords_in = [(x + mx, y + my) for x, y in coords_in]
        coords_out = [(x + mx, y + my) for x, y in coords_out]

        # Create a shapely polygon
        xy = []
        for xy_in, xy_out in zip(coords_in, coords_out):
            xy.append(xy_out)
            xy.append(xy_in)

        return Polygon(xy)


class StarLike(Geometry):
    """Abstract class for star-like geometries."""

    def __init__(self, midpoint: tuple[float, float]) -> None:
        self.midpoint = midpoint

    @abstractmethod
    def radius_at(self, angle: float) -> float:
        """Returns the distance at a given angle to the midpoint."""
        raise NotImplementedError

    def to_polygon(self, refs: int) -> Polygon:
        """Return a shapely polygon representation of self."""
        return self.discretize(refs)

    def discretize(self, refs: int) -> Polygon:
        # TODO: also give angles as alternative.
        angles = [i * 2 * math.pi / refs for i in range(refs)]
        radii = [self.radius_at(angle) for angle in angles]
        mx, my = self.midpoint

        coords = [polar_to_cartesian(angle, rad) for angle, rad in zip(angles, radii)]

        return Polygon([(x + mx, y + my) for x, y in coords])


class Circle(StarLike):
    def __init__(self, midpoint: tuple[float, float], radius) -> None:
        self.radius = radius
        super().__init__(midpoint)

    def radius_at(self, angle: float) -> float:
        return self.radius


class Ellipse(StarLike):
    def __init__(
        self, midpoint: tuple[float, float], axis: tuple[float, float], angle: float
    ) -> None:
        self.axis = axis
        self.angle = angle
        super().__init__(midpoint)

    def radius_at(self, angle: float) -> float:
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
        coefficient: list[tuple[float, float]] = [],
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
        if len(coefficient) == 0:
            coefficient = self.generate_coefficients()

        self.coefficient = coefficient
        self.trig_fun = trigonometric_function(self.coefficient, self.radius)

    def radius_at(self, angle: float) -> float:
        return self.trig_fun(angle)

    def generate_coefficients(self, k: float = 0.2) -> list[tuple[float, float]]:
        """Generate coefficients for the stellar function.

        Args:
            k: Scaling factor.
        """
        j = int(1 / k)
        coeffs_ = [random.uniform(-k, k) for _ in range(j * 2)]

        return list(zip(coeffs_[:j], coeffs_[j:]))


class RainDrop(StarLike):
    """Raindrop geometry."""

    def __init__(self, midpoint: tuple[float, float], a: float, scale: float):
        super().__init__(midpoint)
        if not 0.5 < a and a < 1:
            raise ValueError(f"Value of a must be between 0.5 and 1, but got {a}.")
        self._a = a
        self._scale = scale

    def radius_at(self, angle: float) -> float:
        raise NotImplementedError("Not implemented here")

    def discretize(self, refs):
        angles = [i * 2 * math.pi / refs for i in range(refs)]

        denominators = [
            (1 - self._a + self._a * math.cos(angle)) ** 2
            + self._a**2 * math.sin(angle) ** 2
            for angle in angles
        ]

        x_coords = [
            self._a * math.sin(angle) - self._a * math.sin(angle) / den
            for den, angle in zip(denominators, angles)
        ]
        y_coords = [
            self._a * math.cos(angle) + (1 - self._a + self._a * math.cos(angle)) / den
            for den, angle in zip(denominators, angles)
        ]
        mx, my = self.midpoint
        return Polygon(
            [
                (mx + self._scale * x, my + self._scale * y)
                for x, y in zip(x_coords, y_coords)
            ]
        )


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
    angles = [i * 2 * math.pi / refs for i in range(refs)]

    transformed_coords = rotate_counterclockwise(
        points=[polar_to_cartesian(angle, axis[0], axis[1]) for angle in angles],
        angle=angle,
    )
    mx, my = midpoint

    return Polygon([(x + mx, y + my) for x, y in transformed_coords])


def trigonometric_function(
    coefficients: list[tuple[float, float]], scale: float
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


def rotate_counterclockwise(
    points: list[tuple[float, float]], angle: float
) -> list[tuple[float, float]]:
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
) -> tuple[float, float]:
    """Coordinate transform, from polar to cartesian.

    Args:
        angles: List of angles.
        radii_x: List of radii, for the x-coordinate.
        radii_y: List of radii, for the y-coordinate.
    """
    if rad_y is None:
        rad_y = rad_x
    return (math.cos(angle) * rad_x, math.sin(angle) * rad_y)
