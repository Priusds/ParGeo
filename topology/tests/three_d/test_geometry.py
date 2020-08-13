from topology.three_d.glue import *

def test1():
    # DEFINE GLUE BOX
    lcoords = [
        (0,0,0),
        (2,0,0),
        (2,1,0),
        (0,1,0)
    ]

    depth = 2
    glue = Glue(lcoords, depth)