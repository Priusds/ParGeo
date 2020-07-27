from .star_geometry import *
from .domain import *
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from .utils import rotate_counterclockwise
import copy

def linear(a,b,c,d):
    # returns f:[a,b] -> [c,d], f(x) = mx + n
    # f(a) = c, f(b) = d
    m = (c-d)/(a-b)
    n = c-m*a
    return lambda x: m*x + n


def transform_mesh(mesh, reference_domain, target_domain, rwidth = 0.):
    """
    mesh: is a dolfin mesh

    reference_domain: is a tuple, containing informations about the reference domain
                    reference_domain[0]: rect_ref = (ax_ref, ay_ref)
                    reference_domain[1]: is an object of class Hole

    target_domain: is a tuple, containing informations about the actual domain
                 target_domain[0]: rect_new = (ax_new, ay_new)
                 target_domain[1]: is an object of class Hole

    We assume that the input data is centerd and not rotated. Such that the rect is parallel
    to the origin axis.
    Therefore rect_ref, rect_new only depends on ax, ay. Describing the rectangle [(ax,ay),(-ax,ay),...]

    """

    rect_ref, hole_ref = reference_domain 
    
    angle_new = target_domain.boundary.angle
    midpoint_new = target_domain.boundary.midpoint
    
    hole_new = copy.deepcopy(target_domain.holes[0])
    hole_new.midpoint = (0,0)
    new_mesh = Mesh(mesh)

    rect_new = target_domain.boundary.width/2, target_domain.boundary.height/2


    for ind, p in enumerate(mesh.coordinates()):

        # calculate polar coordinates of p, e.g. angle and l2-distance
        angle = math.atan2(p[1], p[0])
        angle = angle + 2*math.pi if angle < 0 else angle # angle between 0 and 2pi
        p_l2 = np.linalg.norm(p)

        # then we need the l2 distance of the two rects in direction angle
        r_ref = 1/max(abs(np.cos(angle))/rect_ref[0], abs(np.sin(angle))/rect_ref[1])
        r_new = 1/max(abs(np.cos(angle))/rect_new[0], abs(np.sin(angle))/rect_new[1])

        # then we need the l2 distance of the two holes in direction angle
        h_ref = hole_ref.get_radius(angle + angle_new)
        h_ref_agglo = h_ref*(1+rwidth) 
        h_new = hole_new.get_radius(angle + angle_new)
        h_new_agglo = hole_new.get_radius(angle + angle_new)*(1+rwidth)

        #p_l2_hole_dilatation is how much we want to move p according to the dilatation of the hole
        # points on boundary of hole only depend on this
        p_l2_hole_dilatation = linear(h_ref_agglo,r_ref,h_new_agglo,r_new)(p_l2)
        p_l2_hole_dilatation = p_l2_hole_dilatation * np.array([np.cos(angle), np.sin(angle)])

        # p_l2_rect_dilatation is how much we want to move p according to the dilatation of the rectangle
        # points on boundary of rectangle only depend on this
        p_l2_rect_dilatation =  np.array([p[0]*rect_new[0]/rect_ref[0], p[1]*rect_new[1]/rect_ref[1]])

        # Finally for points inbetween hole and rectangle, we choose a convex sum of the factors above
        x = (r_ref - p_l2)/(r_ref - h_ref_agglo)
        
        if x > 1 + 1e-8: # TODO: maybe better error appro
            p_l2_hole_dilatation = linear(h_ref,r_ref,h_new,r_new)(p_l2)
            p_l2_hole_dilatation = p_l2_hole_dilatation * np.array([np.cos(angle), np.sin(angle)])
            x = 1
        p_l2_new = x*p_l2_hole_dilatation + (1-x)*p_l2_rect_dilatation

        # reallocate
        p_l2_new = rotate_counterclockwise(p_l2_new, angle_new)
        
        new_mesh.coordinates()[ind][0] = p_l2_new[0] + midpoint_new[0]
        new_mesh.coordinates()[ind][1] = p_l2_new[1] + midpoint_new[1]

    return new_mesh

def test_transform_mesh():
    """
    Here we test the method above.
    Create some reference domain, its mesh, the target domain and finally using the method above
    create the target mesh.
    """

    lc = 1
    refs = 80
    folder_name = "ref_mesh"
    fname_1 = "S1"
    fname_1 = folder_name + '/' + fname_1
    fname_2 = "S2"
    fname_2 = folder_name + '/' + fname_2

    midpoint = (0,0)
    radius = .8
    circle = Circle(midpoint, radius)
    #stellar = Circle(midpoint, radius + .2)

    stellar = Stellar(midpoint, .5)
    #stellar = Ellipse(midpoint, (.6, .8), 0)

    #stellar.plot()
    #exit()
    #print([stellar.get_radius(x) for x in np.linspace(0,2*np.pi,100)])
    ref_ax = 1
    ref_ay = 1
    ax = 1
    ay = 1


    ref_boundary = [(midpoint[0] - ref_ax, midpoint[1] - ref_ay),
                    (midpoint[0] + ref_ax, midpoint[1] - ref_ay),
                    (midpoint[0] + ref_ax, midpoint[1] + ref_ay),
                    (midpoint[0] - ref_ax, midpoint[1] + ref_ay)]

    ref_topo = Domain(Rectangular(0, 2, 2, (0,0)))
    ref_topo.add_hole(circle)
    ref_topo.writeGeo(fname_1, lc, 10, refs)
    mesh = dolfin_mesh(fname_1)
    ref_domain = ((ref_ax, ref_ay), circle)
    new_domain = ((ax, ay), stellar)


    new_mesh = transform_mesh(mesh, ref_domain, new_domain)
    #exit()
    plt.figure()
    plot(mesh)
    plt.figure()
    plot(new_mesh)
    plt.show()
