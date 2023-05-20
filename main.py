#from bubbles.three_d.bubble import Ellipsoid3D
#import numpy as np
#from bubbles.three_d.sandwich import Sandwich
from bubbles.three_d.glue import Glue
#from bubbles.geo_utils.utils import write_geo
#from bubbles.three_d.bubble import Ellipsoid3D, visualizeLayers, clip
#import pandas as pd
#import os
from bubbles.three_d.bubble import *
from bubbles.geo_utils.utils import write_geo


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
#s1 = Sphere(.1, (1.2,.5,-1))
#glue.add_bubble(s1)



#s3 = StarBubble(.1, (1,.4,-1), None)
#glue.add_bubble(s3,.1)

s4 = StarBubble(.1, (1.8,1,-1.2),None)
glue.add_bubble(s4)

s2 = StarBubble(.1, (.8,1,-1),None)
glue.add_bubble(s2,.5)


# WRITE GEO
glue.to_geo("glue_with_bubbles")







