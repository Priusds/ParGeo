import numpy as np
import matplotlib.pyplot as plt

def rotate_counterclockwise(p, alpha):
    p = np.array(p)
    q = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
    return q @ p


def rotate_clockwise(p, alpha):
    p = np.array(p)
    q_inv = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    return q_inv @ p

def discretize_langloch(langloch, n = 50):
    """
    discretize the langloch
    :param langloch: [m, lenght, width, alpha]
    :return: list of points in ccw
    """
    m, length, width, alpha = langloch
    m = np.array(m)

    #discrete_langloch = []
    r = width/2
    M = int(n/2)
    v = rotate_counterclockwise(np.array([length, 0]), alpha)

    angles_0 = [np.pi * (i/M) + alpha + np.pi/2 for i in range(M + 1)]
    half_circle_0 = [m + r * np.array([np.cos(theta), np.sin(theta)]) for theta in angles_0]
    edge_0 = [half_circle_0[-1] + (i/M)*v for i in range(1, M)]
    angles_1 = [np.pi * (i/M) + alpha + np.pi/2 for i in range(M, 2 * M)]

    new_m = m + v
    half_circle_1 = [new_m + r * np.array([np.cos(theta), np.sin(theta)]) for theta in angles_1]
    edge_1 = [half_circle_1[-1] - (i/M)*v for i in range(1, M)]

    discrete_langloch = half_circle_0 + edge_0 + half_circle_1 + edge_1
    #discrete_langloch = np.array(discrete_langloch)
    #plt.plot(discrete_langloch[:,0], discrete_langloch[:, 1], 'o')
    #plt.show()
    return [(p[0], p[1]) for p in discrete_langloch]

def discretize_ellipse(ellipse, refs):

    """
    discretize the ellipse
    :param ellipse: [m, rads, alpha, lc], n
    :return: list of points in ccw
    """
    m, rads, alpha, lc = ellipse
    n = refs #100

    discrete_ellipse = []

    for i in range(n):
        theta = (i / n) * 2 * np.pi
        discrete_ellipse.append(rotate_counterclockwise(np.array([rads[0] * np.cos(theta), rads[1] * np.sin(theta)]),
                                                        alpha) + m)
    return [(p[0], p[1]) for p in discrete_ellipse]

def stellar_polygon(N, r, center):
    center = np.array(center)
    return [center + r(i/N) * np.array([np.cos((i / N) * 2 * np.pi), np.sin((i / N) * 2 * np.pi)]) for i in range(N)]

def trigonometric_function(coeffients):
    # basis: a_k*cos(k*x) + b_k*sin(k*x)
    """coefficients should be numpy array of dim"""
    n = len(coeffients)
    return lambda x: sum([coeffients[k][0]*np.cos(k*2*np.pi*x) + coeffients[k][1]*np.sin(k*2*np.pi*x) for k in range(n)]) + 4
