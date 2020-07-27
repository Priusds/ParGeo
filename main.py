from topology.three_d.bubble import Ellipsoid3D
import numpy as np
#from topology.three_d.sandwich import Sandwich
from topology.three_d.glue import Glue
from topology.geo_utils.utils import write_geo
from topology.three_d.bubble import Ellipsoid3D, visualizeLayers, clip
import pandas as pd
lcoords = [
    (0,0,0),
    (2,0,0),
    (2,1,0),
    (0,1,0)
]

depth = 2

glue_box = Glue(lcoords, depth)

perm = np.zeros((3,3))
perm[0,0] = 1
perm[1,2] = 1
perm[2,1] = 1


midpoint = (1,1,-1)
radii = (.1,.1,.1)
bubble = Ellipsoid3D(perm @ np.array(midpoint), perm @ np.array(radii))



glue_box.add_bubble(bubble)
glue_box._finish_geometry()

write_geo(glue_box.geo, "test_glue_with_bubble")
exit()
#visualizeLayers(bubble, 40, bubble.getHeightRange(),10)

geo = bubble.get_geo()


def general_rotation(a,b,c):
    R = np.zeros((3,3))
    R[0,0] = np.cos(a)*np.sin(b)
    R[0,1] = np.cos(a)*np.sin(b)*np.sin(c)-np.sin(a)*np.cos(c)
    R[0,2] = np.cos(a)*np.sin(b)*np.cos(c) + np.sin(a)*np.cos(c)

    R[1,0] = np.sin(a)*np.cos(b)
    R[1,1] = np.sin(a)*np.sin(b)*np.sin(c) + np.cos(a)*np.cos(c)
    R[1,2] = np.sin(a)*np.sin(b)*np.cos(c)-np.cos(a)*np.sin(c)

    R[2,0] = -np.sin(b)
    R[2,1] = np.cos(b)*np.sin(c)
    R[2,2] = np.cos(b)*np.cos(c)

    return R


R = general_rotation(.1,.2,.3)
new_points = []
for point in geo["points"]:
    new_point = R @ np.array(point["coords"])
    new_points.append({"id":point["id"], "coords":[c for c in new_point],"lc":point["lc"]})

geo["points"] = new_points

clip_geo = clip(geo, .06)
write_geo(clip_geo, "test_clip")

exit()


write_geo(geo, "bubbly")
write_geo(clip_geo, "clip_bubbly")
exit()

#points_df = pd.DataFrame({{"id":p["id"], "x": p["coords"][0],"y": p["coords"][1],"z":p["coords"][2]} for p in geo["points"]})

#points_df = pd.DataFrame({"id":[p["id"] for p in geo["points"]]})
print("-----------------------------------------------")
points_df = pd.DataFrame(columns=["id","x","y","z","lc"])
for i,point in enumerate(geo["points"]):
    points_df.loc[i]=[point["id"],point["coords"][0],point["coords"][1],point["coords"][2],point["lc"]]
print("-----------------------------------------------")
lines_df = pd.DataFrame(columns=["id","start","end"])
for i,line in enumerate(geo["lines"]):
    lines_df.loc[i] = [line["id"], line["start"], line["end"]]
print("-----------------------------------------------")
lineLoops_df = pd.DataFrame(columns=["id","l1","l2","l3"])
for i,lineLoop in enumerate(geo["lineLoops"]):
    lineLoops_df.loc[i] = [lineLoop["id"], lineLoop["lines"][0], lineLoop["lines"][1], lineLoop["lines"][2]]

# Process the dataFrames
#pd.merge(lineLoops_df, lines_df, on="start")
print("############################################")


print(lineLoops_df.head())
print(lines_df.head())
test = lineLoops_df.set_index("l1").join(lines_df.set_index("id"))
print(test.head())
print(lines_df.loc[lines_df['id'] == 322])
exit()


# Rewrite geo object
new_geo = {"points":[],"lines":[],"lineLoops":[],"surfaces":[]}



write_geo(geo, "bubbly")
write_geo(new_geo, "clip_bubbly")
#bubble.to_geo_old()










exit()

sandwich = Sandwich(1,2,1,5,5)

bubble = Ellipsoid3D((2,2,-2),(.1,.1,.1))

bubble.discretize_full()
bubble.show()
sandwich.add_bubble(bubble)
sandwich.to_geo("sanwich_with_bubble_inside",1)
exit()
midpoint = np.array([0,0,0])
radii = np.array([1,2,3])
angles = np.random.rand(3)*np.pi
E = Ellipsoid3D(midpoint, radii, angles)
h_min, h_max, p_max, p_min = E.getHeightRange()
E.discretize_full()
