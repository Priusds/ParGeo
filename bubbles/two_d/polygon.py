# TODO: Is this deprecated?
import matplotlib.pyplot as plt
import numpy as np

from bubbles.two_d.sat import collide


class Polygon(object):
    def __init__(self, vertices):
        """
        vertices must be in ccw order
        :param vertices:
        """
        self.nv = len(vertices)
        self.vertices = np.array(vertices)
        self.edges = np.array([[vertices[i], vertices[(i + 1) % self.nv]] for i in range(self.nv)])

    def plot(self):
        # just to check if everything is alright
        x = self.vertices[:, 0]
        y = self.vertices[:, 1]
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        plt.plot(x, y)
        plt.show()


def point_on_segment(point, segment):
    """
    point is assumed to be on the line
    """
    p0, p1 = segment
    if p1[0] - p0[0] != 0:
        t = (point[0] - p0[0]) / (p1[0] - p0[0])
    else:
        if p1[1] - p0[1] == 0:
            return p0[0] == point[0] and p0[1] == point[1]
        t = (point[1] - p0[1]) / (p1[1] - p0[1])
    return 0 <= t <= 1


def point_in_polygon(point, polygon):
    point = np.array(point)
    v0, v1 = polygon.edges[0]
    n0 = np.array([v1[1] - v0[1], -(v1[0] - v0[0])])
    sign_ref = np.dot(v0 - point, n0)

    if abs(sign_ref) < 1e-10:
        return not point_on_segment(point, (v0, v1))

    for e in polygon.edges:
        v0, v1 = e
        n0 = np.array([v1[1] - v0[1], -(v1[0] - v0[0])])
        sign = np.dot(v0 - point, n0)

        if sign * sign_ref < 0:
            return False

        if abs(sign) < 1e-10:
            return not point_on_segment(point, (v0, v1))

    return True


def intersection_segment_polygon(seg, poly):
    # convex polygons only
    p0, p1 = seg
    p0 = np.array(p0)
    p1 = np.array(p1)
    t_e = 0
    t_l = 1
    d_s = p1 - p0

    for i, e in enumerate(poly.edges):
        ev = e[1] - e[0]
        n = (ev[1], -ev[0])
        z = - np.dot(p0 - poly.vertices[i], n)
        d = np.dot(d_s, n)

        if abs(d) < 1e-10:

            # then S is parallel to the edge e
            if z < 0:
                # then P_0 is outside of e
                return []
        elif d < 0:
            t = z / d
            # then segment S is entering poly across e
            t_e = max(t_e, t)
            if t_e > t_l:
                return []
        else:
            t = z / d
            t_l = min(t_l, t)
            if t_e > t_l:
                return []

    p_e = p0 + t_e * d_s
    p_l = p0 + t_l * d_s
    output = [p_e, p_l]

    return output


def intersection_poly_poly(polygon0, polygon1):
    return collide(polygon0.vertices, polygon1.vertices)
