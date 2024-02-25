"""Geometry module.

Utility module for creating geometries. The module provides a set of classes
for creating geometries such as rectangles, circles, ellipses, and star-like
geometries. These geometries can be discretized using the `to_polygon` method,
and used in the `Domain` class.

Example:

    >>> from pargeo.domain import Domain
    >>> from pargeo.geometry import Rectangle, Circle
    <BLANKLINE>
    >>> domain = Rectangle(center=(0, 0), width=1, height=1).to_polygon()
    >>> domain = Domain(domain, holes={1})
    >>> hole = Circle(center=(.5, .5), radius=.1).to_polygon()
    >>> domain.add(polygon=hole, level=1)

"""
import math
import random
from abc import ABC, abstractmethod

from pargeo.utils.geometry_utils import (
    discretize_ellipse,
    polar_to_cartesian,
    trigonometric_function,
)
from pargeo.utils.typing_utils import Polygon, Vector


class Geometry(ABC):
    """Abstract class for geometries."""

    @abstractmethod
    def to_polygon(self, *args, **kwargs) -> Polygon:
        """Return a shapely polygon representation of self."""
        pass


class Box(Geometry):
    """Box geometry."""

    def __init__(self, left_bottom: Vector, right_top: Vector):
        """Creates a rectangle.

        Args:
            left_bottom (tuple): The (x, y) coordinates of the left bottom corner.
            right_top (tuple): The (x, y) coordinates of the right top corner.
        """
        self.left_bottom = left_bottom
        self.right_top = right_top
        self.center = tuple(right_top[i] - left_bottom[i] for i in range(2))
        self.width = right_top[0] - left_bottom[0]
        self.height = right_top[1] - left_bottom[1]

    @classmethod
    def from_center(cls, center: Vector, width: float, height: float):
        """
        Creates a box from the center.

        Args:
            center (tuple): The (x, y) coordinates of the center of the box.
            width (float): The width of the box.
            height (float): The height of the box.

        Returns:
            Box: The created box object.
        """
        left_bottom = (center[0] - width / 2, center[1] - height / 2)
        right_top = (center[0] + width / 2, center[1] + height / 2)
        return cls(left_bottom, right_top)

    def to_polygon(self) -> Polygon:
        """Return the corners of the rectangle in counter-clockwise order."""
        p0 = self.left_bottom
        p1 = (self.right_top[0], self.left_bottom[1])
        p2 = self.right_top
        p3 = (self.left_bottom[0], self.right_top[1])

        return Polygon([p0, p1, p2, p3])


class NStar(Geometry):
    """NStar geometry."""

    def __init__(
        self,
        center: Vector,
        radius_in: float,
        radius_out: float,
        n_peaks: int,
        alpha: float = 0.0,
    ):
        """Make an NStar.

        Args:
            center: The center.
            radius_in: The inner radius.
            radius_out: The outer radius.
            N: Number of points.
            alpha: Angle offset. TODO: Explain this.
        """
        if radius_in > radius_out:
            raise ValueError(
                "The parameter setting radius_in > radius_out not allowed."
            )
        self._center = center
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

        # Translate to center
        mx, my = self._center
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

    def __init__(self, center: Vector) -> None:
        self.center = center

    @abstractmethod
    def radius_at(self, angle: float) -> float:
        """Returns the distance at a given angle to the center."""
        raise NotImplementedError

    def to_polygon(self, refs: int) -> Polygon:
        """Return a shapely polygon representation of self."""
        return self.discretize(refs)

    def discretize(self, refs: int) -> Polygon:
        # TODO: also give angles as alternative.
        angles = [i * 2 * math.pi / refs for i in range(refs)]
        radii = [self.radius_at(angle) for angle in angles]
        mx, my = self.center

        coords = [polar_to_cartesian(angle, rad) for angle, rad in zip(angles, radii)]

        return Polygon([(x + mx, y + my) for x, y in coords])


class Circle(StarLike):
    def __init__(self, center: Vector, radius) -> None:
        self.radius = radius
        super().__init__(center)

    def radius_at(self, angle: float) -> float:
        return self.radius


class Ellipse(StarLike):
    def __init__(self, center: Vector, axis: Vector, angle: float = 0) -> None:
        self.axis = axis
        self.angle = angle
        super().__init__(center)

    def radius_at(self, angle: float) -> float:
        # TODO: Check if this is correct
        return math.sqrt(
            (self.axis[0] * math.cos(angle)) ** 2
            + (self.axis[1] * math.sin(angle)) ** 2
        )

    def discretize(self, refs) -> Polygon:
        return discretize_ellipse(self.center, self.axis, self.angle, refs)


class Stellar(StarLike):
    """Stellar bubble."""

    def __init__(
        self,
        center: Vector,
        radius: float,
        coefficient: list[tuple[float, float]] | None = None,
    ):
        """Creates a star-like geometry.

        Args:
            center: The center.
            radius: Approximately the average radius.
            coefficient: numpy array of shape (n, 2) with n > 0.
                Hint: for smooth stellar holes, you can çhoose k from {.1, .2}
                and then `coeffients = np.random.uniform(-k, k, (int(1/k), 2))`
        """
        self.center = center
        self.radius = radius
        if coefficient is None:
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

    def __init__(self, center: Vector, a: float, scale: float):
        super().__init__(center)
        if not (0.5 <= a <= 1):
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
        mx, my = self.center
        return Polygon(
            [
                (mx + self._scale * x, my + self._scale * y)
                for x, y in zip(x_coords, y_coords)
            ]
        )
