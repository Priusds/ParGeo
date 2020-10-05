import numpy as np


def collide(p1, p2):
    """
    p1 and p2 are lists of ordered pairs, the vertices of the polygons in the
    counterclockwise direction.
    :param p1:
    :param p2:
    :return: Return True if the shapes collide. Otherwise, return False.
    """

    p1 = [np.array(v, 'float64') for v in p1]
    p2 = [np.array(v, 'float64') for v in p2]

    edges = edges_of(p1)
    edges += edges_of(p2)

    for o in [orthogonal(e) for e in edges]:
        separates = is_separating_axis(o, p1, p2)
        if separates:
            # they do not collide
            return False

    return True


def centers_displacement(p1, p2):
    """
    Return the displacement between the geometric center of p1 and p2.
    """
    # geometric center
    c1 = np.mean(np.array(p1), axis=0)
    c2 = np.mean(np.array(p2), axis=0)
    return c2 - c1


def edges_of(vertices):
    """
    Return the vectors for the edges of the polygon p.

    p is a polygon.
    """
    edges = []
    n = len(vertices)

    for i in range(n):
        edge = vertices[(i + 1) % n] - vertices[i]
        edges.append(edge)

    return edges


def orthogonal(v):
    """
    Return a 90 degree clockwise rotation of the vector v.
    """
    return np.array([-v[1], v[0]])


def is_separating_axis(o, p1, p2):
    """
    Return True and the push vector if o is a separating axis of p1 and p2.
    Otherwise, return False and None.
    """
    min1, max1 = float('+inf'), float('-inf')
    min2, max2 = float('+inf'), float('-inf')

    for v in p1:
        projection = np.dot(v, o)

        min1 = min(min1, projection)
        max1 = max(max1, projection)

    for v in p2:
        projection = np.dot(v, o)

        min2 = min(min2, projection)
        max2 = max(max2, projection)

    return False if max1 >= min2 and max2 >= min1 else True
