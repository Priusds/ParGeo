from .domain import Domain
from .star_geometry import *
from .minimize import min_problem, min_problem_new
from .logic import *


def cutout(domain, logic):
    """
    Returns the global domain with the cutout boxes and a list of subdomains.
    """
    boundary = domain.boundary
    holes = domain.holes

    bounding_boxes = min_problem_new(boundary, holes) # list of rectangles
    new_domain = Domain(boundary)
    sub_domains = []

    for rect, hole in zip(bounding_boxes, holes):
        # Elements in sub_domains are objects of CRH class

        #new_sub_domain = CRH(rect, hole, nx = 3, ny =3)

        nx, ny = logic.number_boundary_points(rect.width/2, rect.height/2, macro = True)
        sub_rect = PolyRectangular.specify_boundary_nx_ny(rect, nx, ny)

        # Make subdomain
        if hole.midpoint != rect.midpoint:
            fitted_hole = fit_stellar(hole, rect.midpoint, refs = 120, ncoeff = 20)
        new_sub_domain = Domain(sub_rect)
        new_sub_domain.holes.append(fitted_hole)
        sub_domains.append(new_sub_domain)


        new_domain.holes.append(sub_rect)
        domain.bounding_boxes.append(sub_rect)


    return new_domain, sub_domains

def cutout_circle(domain, epsilon):
    # TODO: what is this?
    subdomains = get_subdomains_circle(domain, epsilon)
    new_domain = Domain(domain.boundary)
    subs = []
    for pair in subdomains:
        rect, hole = pair
        new_domain.add_hole(rect)
        subdomain = Domain(rect.boundary)
        subdomain.add_hole(hole)
        subs.append(subdomain)
    return new_domain, subs




def get_subdomains_circle(domain, epsilon):
    # TODO: what is this?
    subdomains = []
    for hole in domain.holes:
        midpoint = hole.midpoint
        dist = hole.radius + epsilon
        p1 = (midpoint[0] - dist, midpoint[1] - dist)
        p2 = (midpoint[0] + dist, midpoint[1] - dist)
        p3 = (midpoint[0] + dist, midpoint[1] + dist)
        p4 = (midpoint[0] - dist, midpoint[1] + dist)

        subdomains.append((Quadrangle([p1, p2, p3, p4]), hole))

    return subdomains
