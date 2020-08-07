#from topology.three_d.bubble import Ellipsoid3D
#import numpy as np
#from topology.three_d.sandwich import Sandwich
from topology.three_d.glue import Glue
#from topology.geo_utils.utils import write_geo
#from topology.three_d.bubble import Ellipsoid3D, visualizeLayers, clip
#import pandas as pd
#import os
from topology.three_d.bubble import *
from topology.geo_utils.utils import write_geo


# DEFINE GLUE BOX
lcoords = [
    (0,0,0),
    (2,0,0),
    (2,1,0),
    (0,1,0)
]

depth = 2
glue = Glue(lcoords, depth)

# DEFINE/ADD BUBBLES
s1 = Sphere(.1, (1.2,.5,-1))
glue.add_bubble(s1,"inside")



s3 = Sphere(.1, (.12,1,-1))
glue.add_bubble(s3,"up")

s4 = Sphere(.1, (1.89,1,-1.2))
glue.add_bubble(s4,"up")

s2 = Sphere(.1, (.8,.5,-1))
glue.add_bubble(s2,"inside")


# WRITE GEO
glue.to_geo("glue_with_bubbles")







