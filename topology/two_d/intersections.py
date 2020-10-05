import numpy as np
from scipy.spatial import ConvexHull


def convex_hull(points):
    new_points = np.array(points)
    hull = ConvexHull(new_points)
    new_points = new_points[hull.vertices]
    return new_points.tolist()


def intersection_hole_corner(hole, edge_0, edge_1, ref):
    val_0, axis_0 = edge_0
    val_1, axis_1 = edge_1

    bool_0, intersections_0 = intersection_hole_edge(hole, edge_0)
    bool_1, intersections_1 = intersection_hole_edge(hole, edge_1)

    if bool_0 and bool_1:
        if ref == "down_left":
            i = intersections_0[0]
            j = intersections_1[0]

            s0 = i if i[0] > val_1 else intersections_0[1]
            s1 = j if j[1] > val_0 else intersections_1[1]
            return True, (s0, s1)

        elif ref == "left_up":
            i = intersections_0[0]
            j = intersections_1[0]

            s0 = i if i[1] < val_1 else intersections_0[1]
            s1 = j if j[0] > val_0 else intersections_1[1]
            return True, (s0, s1)

        elif ref == "up_right":
            i = intersections_0[0]
            j = intersections_1[0]

            s0 = i if i[0] < val_1 else intersections_0[1]
            s1 = j if j[1] < val_0 else intersections_1[1]
            return True, (s0, s1)

        else:
            i = intersections_0[0]
            j = intersections_1[0]

            s0 = i if i[1] > val_1 else intersections_0[1]
            s1 = j if j[0] < val_0 else intersections_1[1]
            return True, (s0, s1)

    else:
        return False, None


def intersection_hole_edge(hole, edge):
    hole_shifted = hole[1:] + [hole[0]]
    segments = zip(hole, hole_shifted)

    counter = 0
    intersections = []
    for seg in segments:
        bool_, s_points = intersection_segment_edge(seg, edge)

        if bool_:
            counter += 1
            if counter > 2:
                return False, None
            intersections.append(s_points)

    if len(intersections) != 2:
        return False, None
    return True, intersections


def intersection_segment_edge(segment, edge):
    # here we calculate the intersection points
    val, axis = edge
    p = segment[0]
    q = segment[1]

    if axis == 'x':
        bool1 = p[0] < val < q[0]
        bool2 = p[0] > val > q[0]

        if bool1 or bool2:
            # y = mx + n
            assert q[0] - p[0] != 0
            m = (q[1] - p[1]) / (q[0] - p[0])
            n = p[1] - m * p[0]
            # intersection = (p[1] + q[1]) / 2
            intersection = m * val + n
            return True, (val, intersection)
        return False, None
    else:
        bool1 = p[1] < val < q[1]
        bool2 = p[1] > val > q[1]

        if bool1 or bool2:
            # y = mx + n
            assert q[0] - p[0] != 0
            m = (q[1] - p[1]) / (q[0] - p[0])
            n = p[1] - m * p[0]
            # intersection = (p[0] + q[0]) / 2
            intersection = (val - n) / m
            return True, (intersection, val)
        return False, None
