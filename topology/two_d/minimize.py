from .star_geometry import *
import numpy as np
from scipy.optimize import minimize as scipy_minimize
from .domain import Domain
from shapely.geometry import LinearRing as lr


dist_hole_hole = lambda hole0, hole1, refs0=50, refs1=50: lr(hole0.discretize_hole(refs0)).distance(lr(hole1.discretize_hole(refs1)))
dist_hole_poly = lambda hole, list_, refs=50: lr(hole.discretize_hole(refs)).distance(lr(list_))

def min_problem(boundary, holes):
    boundary_points = boundary.discretize_hole(refs = 50)
    n = len(holes)
    #holes = map(lambda x: np.array(x), holes)

    def f(angles):

        rects = list(map(lambda X: get_midpoint_ax_ay(X[0], X[1]), zip(holes,angles)))
        res = 0.
        dists_E = [dist_hole_poly(r, boundary_points)**2 for r in rects]
        dists_R = [dist_hole_hole(rects[i], rects[j])**2 for i in range(len(rects)) for j in range(i + 1, len(rects))]

        if 0. in  dists_E:
            res = np.inf
        else:
            res += sum(1. / np.array(dists_E))

        if 0. in dists_R:
            res = np.inf
        else:
            res += sum(1. / np.array(dists_R))


        #print("================================")
        #print("f({a}) = {b}".format(a=angles,b=res))
        #print([round(de, 2) for de in dists_E])
        #print([round(dr,2) for dr in dists_R])
        #print("minE = ", min(dists_E))
        #print("minR = ", min(dists_R)) if len(dists_R) > 0 else None
        #print("================================")
        return res

    optimize = False


    if optimize:
        #res = scipy_minimize(f, np.zeros(n), method = 'SLSQP', constraints = {'type':'ineq', 'fun':cons0})

        #N = 100

        #for _ in range(N):
        #    alphas = np.pi * np.random.rand(n)
        #    if f(alphas) < np.inf:
        #        break


        x0 = np.pi * np.random.rand(n)

        res = scipy_minimize(f, x0, method = 'Powell', options = {'maxiter':10, 'maxfev' : 10, 'disp' : True})
        optimal_angles = res.x

    else:
        N = 100
        alpha_opt = np.zeros(n)
        best = f(alpha_opt)
        for _ in range(N):
            alpha = np.pi*np.random.rand(n)/2

            new_f = f(alpha)

            if new_f < best:
                alpha_opt = alpha
                best = new_f

        optimal_angles = alpha_opt

    

    rects = list(map(lambda X: get_midpoint_ax_ay(X[0], X[1]), list(zip(holes, optimal_angles))))
    if not valid(rects, boundary_points):
        optimal_angles = np.zeros(n)
        rects = list(map(lambda X: get_midpoint_ax_ay(X[0], X[1]), list(zip(holes, optimal_angles))))
        #raise ValueError("Optimization problem did't work")
    print("++++++++RESULT+++++++++++++")
    f(optimal_angles)

    # make sure angle to be between 0 and pi/2, and ax >= ay
    if min(optimal_angles) > np.pi/2:
        raise ValueError("optimal angle {0} should be between 0 and pi/2".format(min(optimal_angles)))

    for rect in rects:
        q, ax, ay, m = rect.angle, rect.width/2, rect.height/2, rect.midpoint
        if ax < ay:
            temp = ax
            ax = ay
            ay = temp
            q = q + np.pi/2
            rect = Rectangular(angle = q, width = ax*2, height = ay*2, midpoint = m)

    return rects



def min_problem_new(boundary, holes):
    refs = 50
    discrete_holes = [hole.discretize_hole(refs) for hole in holes]
    nholes = len(holes)
    angles = np.zeros(nholes)
    
    niter = 0
    while niter < 1000:
        niter += 1
        rects = list(map(lambda X: get_rect_new(X[0], X[1]), list(zip(discrete_holes, angles))))
        collisions = []
        for i in range(nholes):
            if dist_hole_hole(boundary, rects[i]) <= 1e-6:
                collisions.append([i])
            for j in range(i+1, nholes):
                if dist_hole_hole(rects[i], rects[j]) <= 1e-6:
                    collisions.append((i,j))
        if len(collisions) == 0:
            return rects

        vertices = set([v for c in collisions for v in c])
        for i in vertices:
            angles[i] = np.random.uniform(0, np.pi/2)
    
    #return list(map(lambda X: get_rect_new(X[0], X[1]), list(zip(discrete_holes, angles))))
    raise ValueError("Optimization didn't work")


def valid(rects, boundary_points):
    """
    Checks if this boundary with rects is a valid domain.
    """
    minimal_boundary_dist = min([dist_hole_poly(r, boundary_points) for r in rects]) if len(rects) > 0 else Domain.NU
    minimal_rect_to_rect_dist = min([dist_hole_hole(rects[i], rects[j]) for i in range(len(rects)) for j in range(i + 1, len(rects))]) if len(rects) > 1 else Domain.GAMMA

    if minimal_boundary_dist >= Domain.NU and minimal_rect_to_rect_dist >= Domain.GAMMA:
        return True
    return False


def get_rect_new(hole_discrete, alpha):
    if not 0 <= alpha <= np.pi/2:
        raise ValueError("Optimal angle {0} of cutout rectangle should be between 0 and pi/2".format(alpha))
    epsilon = Domain.EPSILON

    
    points = np.array(hole_discrete)

    w1 = np.array([np.cos(alpha), np.sin(alpha)])
    w2 = np.array([np.cos(alpha + np.pi/2), np.sin(alpha + np.pi/2)])
    w3 = np.array([-np.cos(alpha), -np.sin(alpha)])
    w4 = np.array([-np.cos(alpha + np.pi/2), -np.sin(alpha + np.pi/2)])
    wl = np.array([w1,w2,w3,w4]).transpose()

    scalar_prods = points @ wl
    extremal_points_index = np.argmax(scalar_prods, axis=0)
    max_scalar_prods = [scalar_prods[:,i][extremal_points_index[i]] for i in range(4)]
    width = max_scalar_prods[0] + max_scalar_prods[2]
    height = max_scalar_prods[1] + max_scalar_prods[3]
    midpoint = w1*(max_scalar_prods[0] - max_scalar_prods[2])/2 + w2*(max_scalar_prods[1] - max_scalar_prods[3])/2

    ax = 0.5*width*(1+2*epsilon)
    ay = 0.5*height*(1+2*epsilon)
    midpoint = midpoint.tolist()

    if ax < ay:
        temp = ax
        ax = ay
        ay = temp
        alpha = alpha + np.pi/2

    return Rectangular(angle = alpha, width = ax*2, height = ay*2, midpoint = midpoint)


def get_midpoint_ax_ay(hole, alpha):
    if not 0 <= alpha <= np.pi/2:
        raise ValueError("Optimal angle {0} of cutout rectangle should be between 0 and pi/2".format(alpha))
    epsilon = Domain.EPSILON

    refs = 50
    points = np.array(hole.discretize_hole(refs))

    w1 = np.array([np.cos(alpha), np.sin(alpha)])
    w2 = np.array([np.cos(alpha + np.pi/2), np.sin(alpha + np.pi/2)])
    w3 = np.array([-np.cos(alpha), -np.sin(alpha)])
    w4 = np.array([-np.cos(alpha + np.pi/2), -np.sin(alpha + np.pi/2)])
    wl = np.array([w1,w2,w3,w4]).transpose()

    scalar_prods = points @ wl
    extremal_points_index = np.argmax(scalar_prods, axis=0)
    max_scalar_prods = [scalar_prods[:,i][extremal_points_index[i]] for i in range(4)]
    width = max_scalar_prods[0] + max_scalar_prods[2]
    height = max_scalar_prods[1] + max_scalar_prods[3]
    midpoint = w1*(max_scalar_prods[0] - max_scalar_prods[2])/2 + w2*(max_scalar_prods[1] - max_scalar_prods[3])/2

    ax = 0.5*width*(1+2*epsilon)
    ay = 0.5*height*(1+2*epsilon)
    midpoint = midpoint.tolist()

    if ax < ay:
        temp = ax
        ax = ay
        ay = temp
        alpha = alpha + np.pi/2

    return Rectangular(angle = alpha, width = ax*2, height = ay*2, midpoint = midpoint)
