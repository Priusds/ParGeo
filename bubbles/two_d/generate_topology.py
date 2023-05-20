# TODO: Check if this module is deprecated.
from .domainWithHoles import DomainWithHoles as dwh
from .hole import Circle
import json
import random as rd


def generate_circle_topology(file_name, n_holes, radius, width = 9, height = 9, dir = None):
    rect_out = [(0, 0), (width, height)]
    rect_in = [(width/3, height/3), (2*width/3, 2*height/3)]
    topology = dwh(file_name, rect_out, 1, rect_in, 1, dir = dir) # the 1 for lc_in and lc_out doesn't play a role

    i = 0
    while i < n_holes:
        midpoint_x = rd.uniform(0, width)
        midpoint_y = rd.uniform(0, height)
        midpoint = (midpoint_x, midpoint_y)

        circle = Circle({"type":"circle", "radius":radius,"midpoint" : midpoint})
        if topology.new_add_hole(circle):
            i = i + 1

    with open(file_name + ".json", 'w') as fp:
        json.dump(topology.valid_list, fp, indent=1)


        #with open('data.json', 'r') as fp:
        #data = json.load(fp)
