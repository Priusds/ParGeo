from abc import ABC, abstractmethod

import matplotlib.pyplot as plt
import numpy as np


class Hole(ABC):
    def __init__(self, type_):
        self.type = type_

    @abstractmethod
    def discretize_hole(self):
        """returns a list of points (vertices) in ccw order"""
        raise NotImplementedError

    @abstractmethod
    def to_dict(self):
        """returns the hole informations in dict format"""
        raise NotImplementedError

    def info(self):
        """prints informations about the hole"""
        print(self.to_dict())

    def getPoint(self, angle):
        """returns the corresponding point on the submanifold for given angle"""
        raise NotImplementedError


class Circle(Hole):
    def __init__(self, midpoint, radius):
        super().__init__("circle")
        self.midpoint = midpoint
        self.radius = radius

    def to_dict(self):
        return {"type": self.type, "midpoint": self.midpoint, "radius": self.radius}

    def discretize_hole(self, refs):
        return discretize_ellipse([self.midpoint, (self.radius, self.radius), 0], refs)

    def getPoint(self, angle, axis):
        if axis == 0:
            return self.radius * np.cos(angle) + self.midpoint[0]
        elif axis == 1:
            return self.radius * np.sin(angle) + self.midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Ellipse(Hole):
    def __init__(self, midpoint, axis, angle):
        super().__init__("ellipse")
        self.midpoint = midpoint
        self.axis = axis
        self.angle = angle

    def to_dict(self):
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "axis": self.axis,
            "angle": self.angle,
        }

    def discretize_hole(self, refs):
        return discretize_ellipse([self.midpoint, self.axis, self.angle], refs)

    def getPoint(self, angle, axis):
        if axis == 0:
            return self.axis[0] * np.cos(angle) + self.midpoint[0]
        elif axis == 1:
            return self.axis[1] * np.sin(angle) + self.midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Stellar(Hole):
    """
    coeffient must be of type numpy.ndarray
    For smooth forms:
    coeffients = np.random.uniform(-k, k, (j, 2))
    where k in [.1, .2] and j = 1 / k

    """

    def __init__(self, midpoint, radius, coefficient=None):
        super().__init__("stellar")
        self.midpoint = midpoint
        self.coefficient = (
            np.random.uniform(-0.2, 0.2, (2, 2)) if coefficient is None else coefficient
        )
        self.radius = radius

        self.f = trigonometric_function(self.coefficients, self.radius)

    def to_dict(self):
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "coefficient": self.coefficient.tolist(),
            "radius": self.radius,
        }

    def discretize_hole(self, refs):
        return discetize_stellar_polygon(
            self.midpoint, self.radius, self.coefficient, refs
        )

    def getPoint(self, angle, axis):
        r = self.f(angle)

        if axis == 0:
            return r * np.cos(angle) + midpoint[0]
        elif axis == 1:
            return r * np.sin(angle) + midpoint[1]
        else:
            raise ValueError("Only axis = 0 or 1 supported")


class Rectangular(Hole):
    def __init__(self, midpoint, width, height):
        super().__init__("rectangular")
        self.midpoint = midpoint
        self.width = width
        self.height = height

    def to_dict(self):
        return {
            "type": self.type,
            "midpoint": self.midpoint,
            "width": self.width,
            "height": self.height,
        }

    def discretize_hole(self, refs):
        midpoint = self.midpoint
        width = self.width
        height = self.height
        p0 = (midpoint[0] - width / 2, midpoint[1] - height / 2)
        p1 = (midpoint[0] + width / 2, midpoint[1] - height / 2)
        p2 = (midpoint[0] + width / 2, midpoint[1] + height / 2)
        p3 = (midpoint[0] - width / 2, midpoint[1] + height / 2)

        return [p0, p1, p2, p3]


def discretize_ellipse(ellipse, refs):
    """
    discretize the ellipse
    :param ellipse: [m, rads, alpha], n
    :return: list of points in ccw
    """
    m, rads, alpha = ellipse
    n = refs  # 100
    discrete_ellipse = []

    for i in range(n):
        theta = (i / n) * 2 * np.pi
        discrete_ellipse.append(
            rotate_counterclockwise(
                np.array([rads[0] * np.cos(theta), rads[1] * np.sin(theta)]), alpha
            )
            + m
        )
    return [(p[0], p[1]) for p in discrete_ellipse]


def discretize_ellipse_noisy(ellipse, refs, var):
    """
    var ~ radius/10
    discretize the ellipse with additional noise
    :param ellipse: [m, rads, alpha], n
    :return: list of points in ccw
    """
    m, rads, alpha = ellipse
    n = refs  # 100
    discrete_ellipse = []
    noise = np.random.normal(0, var, refs)

    for i in range(n):
        theta = (i / n) * 2 * np.pi
        discrete_ellipse.append(
            rotate_counterclockwise(
                np.array(
                    [
                        (rads[0] + noise[i]) * np.cos(theta),
                        (rads[1] + noise[i]) * np.sin(theta),
                    ]
                ),
                alpha,
            )
            + m
        )
    return [(p[0], p[1]) for p in discrete_ellipse]


def discetize_stellar_polygon(midpoint, radius, coefficients, refs):
    X = np.linspace(0, 2 * np.pi, refs, endpoint=False)
    f = trigonometric_function(coefficients, radius)
    rad = [f(x) for x in X]

    polygon = [
        (r * np.cos(omega) + midpoint[0], r * np.sin(omega) + midpoint[1])
        for r, omega in zip(rad, X)
    ]

    return polygon


def trigonometric_function(coefficients, radius):
    # basis: a_k*cos(k*x) + b_k*sin(k*x)
    """coefficients should be numpy array of shape n x 2"""

    n = len(coefficients)
    return lambda x: radius * np.exp(
        sum(coefficients[0]) / 2
        + sum(
            [
                coefficients[k][0] * np.cos(k * x) + coefficients[k][1] * np.sin(k * x)
                for k in range(1, n)
            ]
        )
    )


def rotate_counterclockwise(p, alpha):
    p = np.array(p)
    q = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
    return q @ p


def rotate_clockwise(p, alpha):
    p = np.array(p)
    q_inv = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    return q_inv @ p
