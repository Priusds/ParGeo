import numpy as np
from scipy.optimize import LinearConstraint, minimize, root_scalar
from scipy.spatial import ConvexHull


def convex_hull(points):
    new_points = np.array(points)
    hull = ConvexHull(new_points)
    new_points = new_points[hull.vertices]
    return new_points.tolist()


def intersection_hole_corner(hole, edge_0, edge_1, ref, method="linear"):
    val_0, axis_0 = edge_0
    val_1, axis_1 = edge_1

    bool_0, intersections_0 = intersection_hole_edge(hole, edge_0, method)
    bool_1, intersections_1 = intersection_hole_edge(hole, edge_1, method)

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


def intersection_hole_edge(hole, edge, method="linear"):
    """
    hole is a tuple of (dhole, chole)
        dhole is discretization of the underlying continuous hole chole.
    """
    dhole = hole[0]
    # print("dhole = ", dhole)

    if method == "linear":
        dhole_shifted = dhole[1:] + [dhole[0]]
        segments = zip(dhole, dhole_shifted)

        counter = 0
        intersections = []
        for seg in segments:
            # print("seg : ", seg)
            bool_, s_points = intersection_segment_edge(seg, edge)

            if bool_:
                counter += 1
                if counter > 2:
                    return False, None
                intersections.append(s_points)

        if len(intersections) != 2:
            return False, None
        return True, intersections

    elif method == "optimise":
        # TODO if stellar is special form like circle / ellipse make analytical intersections

        if not edge[1] in ["x", "y"]:
            raise ValueError("Non proper edge input.")

        val, axis = edge[0], 0 if edge[1] == "x" else 1

        chole = hole[1]

        s = 2 * np.pi / len(dhole)
        thetas = [i * s for i in range(len(dhole))]
        thetas_shifted = thetas[1:] + [thetas[0]]

        # dhole = hole[0]
        # dhole_shifted = dhole[1:] + [dhole[0]]
        # pointData = zip(dhole,dhole_shifted, theta, theta_shifted)
        thetaPairs = zip(thetas, thetas_shifted)

        def f(t):
            # print("a = {a}, b = {b}".format(a= chole.getPoint(t, axis ), b = val))
            return chole.getPoint(t, axis) - val

        counter = 0
        intersections = []
        for tp in thetaPairs:
            lb, ub = tp

            if (
                f(lb) * f(ub) <= 0
            ):  # bool_, s_points = intersection_segment_edge(seg, edge)
                res = root_scalar(f, bracket=[lb, ub])
                counter += 1
                P = [chole.getPoint(res.root, 0), chole.getPoint(res.root, 1)]
                # print("projection error : ", abs(P[axis] - val))

                # print("P = {P}".format(P=P))

                P[
                    axis
                ] = val  # project the axis coordinate to the correct retangular edge height

                # print("P projected = {P}".format(P=P))
                intersections.append(P)

        # print("LEN P = ", len(intersections))
        if len(intersections) != 2:
            return False, None
        return True, intersections

    else:
        raise NotImplementedError("Other methods not avaiable.")


def intersection_segment_edge(segment, edge):
    # here we calculate the intersection points
    val, axis = edge
    p = segment[0]
    q = segment[1]

    if axis == "x":
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
