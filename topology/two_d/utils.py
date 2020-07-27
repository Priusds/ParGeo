import numpy as np
import math
cos = math.cos
sin = math.sin

# TODO: Remove old code and add everything related to geometry

def angle(v):
    """Given a 2d-vector returns angle with (1,0) vector."""
    # TODO: Check correctness
    angle = math.atan2(v[1],v[0]) 
    return angle if angle > 0 else 2*math.pi + angle

def correct(new_angle):
    # make sure 0 <= new_angle < 2PI
    new_angle = new_angle if new_angle >= 0 else 2*np.pi + new_angle 
    new_angle = new_angle if new_angle < 2*np.pi else new_angle - 2*np.pi*int(new_angle/(2*np.pi))
    return new_angle
    
def get_rect(z):
    """Returns list of boundary points, for given data: (angle, ax, ay, midpoint)"""
    phi, a_x, b_y, midpoint = z

    midpoint = np.array(midpoint)
    
    rect = np.array([[-1, 1, 1, -1],
                    [-1, -1, 1, 1]])

    d = np.diag([a_x, b_y])

    q = rot_matrix(phi)

    return ((q @ d @ rect).transpose() + midpoint).tolist()


def rotate_counterclockwise(p, alpha, midpoint = (0,0)):
    cp = p[0] - midpoint[0], p[1] - midpoint[1]
    return cp[0]*cos(alpha)-cp[1]*sin(alpha), cp[0]*sin(alpha) + cp[1]*cos(alpha)

def rotate_clockwise(p, alpha, midpoint = (0,0)):
    cp = p[0] - midpoint[0], p[1] - midpoint[1]
    return cp[0]*cos(alpha) + cp[1]*sin(alpha), -cp[0]*sin(alpha) + cp[1]*cos(alpha)

def rot_matrix(alpha):
    # counterclockwise
    return np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])

def get_rect_from_boundary(boundary):
    """
    returns midpoint, angle, ax, ay of rectangle given it's vertices

    ax >= ay always.

    """
    boundary_np = np.array(boundary)
    midpoint_np = np.mean(boundary_np, axis = 0)

    ax = np.linalg.norm(boundary_np[0] - boundary_np[1])/2
    ay = np.linalg.norm(boundary_np[1] - boundary_np[2])/2

    x = boundary_np[1] - boundary_np[0]
    angle = math.atan2(x[1], x[0])
    angle = angle + 2*math.pi if angle < 0 else angle # angle must be positive
    angle = angle if angle <= math.pi else angle - math.pi

    if ax < ay:
        temp = ay
        ay = ax
        ax = temp
        angle = angle + math.pi/4 if angle <= math.pi/4 else angle - math.pi/4

    return midpoint_np.tolist(), angle, ax, ay


