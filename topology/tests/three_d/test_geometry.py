from topology.three_d.glue import *
from topology.three_d.clip import *
from topology.three_d.utils import locate
from topology.three_d.bubble import Sphere,clip
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

def test2():

    n = (0,1,0)
    s = (66,1,12)
    p = (4,5,6)


    e = EasyClip(EuclideanProjection(n,s))
    print(e.projectPoint(p))

def test3():
    s = Sphere(1,(0,0,0))
    accuracy = .1
    geo = s.discretize(accuracy)
    n = (-10,0,0)
    s = (0.5,0,0)
    clippingPlane = (n,s)


    #print("TEYT", locate(point, n,s)) 
    #print(locate(npoint, n, s))
    #print(sum([(npoint[i] - s[i])*n[i] for i in range(3)]))
    newgeo = clip(geo, clippingPlane, EuclideanProjection, EasyClip)
   # print(newgeo["points"])
    write_geo(newgeo, "tittifucka")

def test4():

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
    s1 = Sphere(.1, (1,.05,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (1,.95,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (1.95,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (.01,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (1,.5,-1.95)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (1,.5,-.01)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)

    s1 = Sphere(.1, (1,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    glue.add_bubble(s1)
   # s1 = Sphere(.1, (1,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    #s1 = Sphere(.1, (.5,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)


    #s3 = StarBubble(.1, (1,.4,-1), None)
    #glue.add_bubble(s3,.1)

    #s4 = StarBubble(.1, (1.8,1,-1.2),None)
    #glue.add_bubble(s4)

    #s2 = StarBubble(.1, (.8,1,-1),None)
    #glue.add_bubble(s2,.5)


    # WRITE GEO
    glue.to_geo("glue_with_bubbles")




