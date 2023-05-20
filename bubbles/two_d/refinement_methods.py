# Methods for setting the hole refinement, for now only depending on the distance from center
# Refinement methods should have following format:
# Input: rect_out, rect_in, hole_infos, decay
# rect_out = [(0,0), (width, height)]
# rect_in = [(width/3, height/3), (2*width/3, 2*height/3)]
# hole_infos is hole specific, for example midpoint for circle
# decay is a parameter for tuning the decay
import numpy as np


def linear_decay(rect_out, rect_in, midpoint, decay=20):
    # extract infos
    midpoint_x, midpoint_y = midpoint
    width = rect_out[1][0]
    height = rect_out[1][1]

    cdom = (width / 2, height / 2)  # center of the domain, used for refinement level
    rad_subarea = height / 6

    # distance from subarea
    dist_cm = max(abs(cdom[0] - midpoint_x), abs(cdom[1] - midpoint_y))

    nrefs = 50 if dist_cm <= rad_subarea else 50 - decay * (dist_cm - rad_subarea)

    return round(nrefs)


def sigmoid_decay(rect_out, rect_in, midpoint, decay=8):
    """sigmoid dacay for rect_in quadratic centerd in rect_out quadratic"""
    # extract infos
    midpoint_x, midpoint_y = midpoint
    width_out = rect_out[1][0]
    width_in = rect_in[1][0] - rect_out[0][0]
    inner_rad = width_in / 2

    domain_center = (
        width_out / 2,
        width_out / 2,
    )  # center of the domain, used for refinement level

    # distance from subarea
    dist = max(abs(domain_center[0] - midpoint_x), abs(domain_center[1] - midpoint_y))
    # scale distance, starting on the boundary
    dist = dist - inner_rad if dist > inner_rad else 0

    # max distance
    max_dist = domain_center[0] - inner_rad

    assert 0 <= dist / max_dist <= 1
    # refeinement of the hole will be between 50 and 50 - max_ref (Default = 20)
    max_ref = 30
    nrefs = 50 - sigmoid(dist / max_dist, max_ref, decay)
    assert nrefs > 3
    return int(round(nrefs))


def sigmoid(x, max_val, alpha):
    """x shoud be between 0 and 1, the output between 0 and max"""
    x = alpha * x
    y = 2 * (np.exp(x) / (1 + np.exp(x)) - 0.5)
    return max_val * y
