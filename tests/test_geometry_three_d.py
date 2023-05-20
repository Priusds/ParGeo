from bubbles.three_d.glue import *
from bubbles.three_d.clip import *
from bubbles.three_d.utils import locate, merge_loops
from bubbles.three_d.bubble import Sphere,clip, Ellipse
from bubbles.three_d.utils import locate
from bubbles.three_d.bubble import Sphere,clip, _targetPoint
#from topology.three_d.geometry import *
from bubbles.three_d.layers import *
from bubbles.three_d.merge import   Sandwich, makeSomeSandwiches, mergeGeometries, Blade, makeSomeSandwiches2


import mystic as my
# tools
from mystic.monitors import VerboseMonitor
from mystic.math.measures import mean, impose_mean
from mystic.math import almostEqual
from mystic.penalty import quadratic_equality

import numpy as np


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
  #  s1 = Sphere(.1, (1,.05,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
  #  glue.add_bubble(s1)

   # s1 = Sphere(.1, (1,.95,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
   # glue.add_bubble(s1)

    #s1 = Sphere(.1, (1.95,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
   # glue.add_bubble(s1)

    #s1 = Sphere(.1, (.01,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    #s1 = Sphere(.1, (1,.5,-1.96)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    #s1 = Sphere(.1, (.5,.5,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    #s1 = Sphere(.1, (1,1.97,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1,1)

   # s1 = Sphere(.1, (0.01,.05,-.5)) # TODO: check problem with Sphere(.3, (1,"1",-1))
   # glue.add_bubble(s1)

   # s1 = Sphere(.1, (0.01,.05,-1.5)) # TODO: check problem with Sphere(.3, (1,"1",-1))
   # glue.add_bubble(s1)


    #s1 = Sphere(.1, (.5,.05,-1)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    #s1 = Sphere(.1, (1.2,.05,-.2)) # TODO: check problem with Sphere(.3, (1,"1",-1))
    #glue.add_bubble(s1)

    

    #s3 = StarBubble(.1, (1,.4,-1), None)
    #glue.add_bubble(s3,.1)

    #s4 = StarBubble(.1, (1.8,1,-1.2),None)
    #glue.add_bubble(s4)

    #s2 = StarBubble(.1, (.8,1,-1),None)
    #glue.add_bubble(s2,.5)


    # WRITE GEO
    glue.to_geo("glue_with_bubbles")




def test5():
    s1 = Sphere(1, (0,0,0))
    geoB = s1.discretize(.3)
    h1 = (0,0,1)
    h2 = (1,1,0)
    h3 = (0,1,0)
    geo0, loop0 = clip(geoB, (h1,(0,0,0)), pointProjection = EuclideanProjection,
                            clippingRule = EasyClip)
   
    geo1, loop1 = clip(geoB, (h2,(0,0,0)), pointProjection = EuclideanProjection,
                            clippingRule = EasyClip)

    p0 = [geoB["lines"][abs(l)][0] for l in loop0[1]]
    p1 = [geoB["lines"][abs(l)][0] for l in loop1[1]]
    print(set(p0).intersection(set(p1)))
    #write_geo(geo, "test")

def test6():
    s1 = Sphere(1, (0.05,0.05,0.05))
    geoB = s1.discretize(.7)
    h0 = (0,0,1.6)
    h1 = (1.6,0,0)

    geo0,loop0 = clip(geoB, (h0,(0,0,0)), pointProjection = EuclideanProjection,
                            clippingRule = EasyClip)
    #print(loop0)
    geo1, loop1 = clip(geo0, (h1,(0.1,0,0)), pointProjection = EuclideanProjection,
                            clippingRule = EasyClip)
   # print(loop1)



    p0 = [geoB["lines"][abs(l)][0] for l in loop0[1]]+[geoB["lines"][abs(l)][1] for l in loop0[1]]
    p1 = [geoB["lines"][abs(l)][0] for l in loop1[1]]+[geoB["lines"][abs(l)][1] for l in loop1[1]]
    #print(set(p0).intersection(set(p1)))

    loop0_clip0 = []
    loop0_clip1 = []
    out = False
    for i,l in enumerate(loop0[1]):
        if abs(l) in geo1["lines"] and not out:
            loop0_clip0.append(l)
        elif abs(l) in geo1["lines"] and out:
            loop0_clip1.append(l)
        else:
            out = True
    loop0 = loop0_clip1 + loop0_clip0

    print(loop0)
    write_geo(geo1, "test")




def testProjectionClipping():

    s1 = Sphere(1, (0,0,0))
    geo = s1.discretize(0.1)

    myProjection= NumericPointProjection
    myClipping = EasyClip
    clippingPlane = ((0,1,0), (0,0.5,0))

    phi = lambda p : np.array(_targetPoint(p[0],p[1],s1.trafo(p[0],p[1])))

    geo2, loop = clip(geo, clippingPlane   = clippingPlane,
                          pointProjection = myProjection,
                          clippingRule    = myClipping,
                          args = phi)
    write_geo(geo2, "testClipping")



    myProjection= EuclideanProjection
    myClipping = EasyClip
    clippingPlane = ((0,1,0), (0,0.5,0))

    geo3, loop = clip(geo, clippingPlane   = clippingPlane,
                          pointProjection = myProjection,
                          clippingRule    = myClipping,
                          args = None)
    write_geo(geo3, "testClipping_eucl")



                                                  


def test8():
    w = Waffle(lcoords = [(0,0,0),(2,0,0)],
                depth = 3,
                nLayers = 5,
                hightLayer = .5)
    w.to_geo("Waffle")





def test11():
    # Make Glue
    depth = 2
    x = 0
    y = 0
    z = 0

    width = 2
    hight = 1

    lcoords = [
        (x,y,z),
        (x,y-width,z),
        (x+hight,y-width,z),
        (x+hight,y,z)
    ]

    

    glue = Glue(lcoords, depth)
    glue._finish_geometry()
    write_geo(glue.geometry, "titten")

def test12():
    # Rotate sandwich
    depth = 2
    x = 0
    y = 0
    z = 0

    width = 2
    hight = 1

    lcoords = [
        (x,y,z),
        (x+width,y,z),
        (x+width,y+hight,z),
        (x,y+hight,z)
    ]

    

    glue = Glue(lcoords, depth)



    glue._finish_geometry()

    w1 = Waffle(lcoords = [(x,y-3*.1,z),(x+width,y-3*.1,z)],
                depth = depth,
                nLayers = 3,
                hightLayer = .1)

    w2 = Waffle(lcoords = [(x,y+hight,z),(x+width,y+hight,z)],
                depth = depth,
                nLayers = 1,
                hightLayer = .1)
    
    sandwich = Sandwich(w1, glue, w2).geometry

    for pointID, pointDict in sandwich["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (-pointCoords[1], pointCoords[0], pointCoords[2])
        sandwich["points"][pointID]["coords"] = newCoords
        #print(sandwich["points"][pointID])
    write_geo(sandwich, "julubu")

def test13():
    # Make Layer
    depth = 2
    x = 0
    y = 0
    z = 0

    width = 2
    hight = 1
    w1 = Waffle(lcoords = [(x,y-3*.1,z),(x+width,y-3*.1,z)],
                depth = depth,
                nLayers = 3,
                hightLayer = .1)
    write_geo(w1.geometry, "waffleWonder")



def test16():
    depth = 2
    x = 1
    y = 1
    z = 0

    width = 2
    hight = 1

    lcoords = [
        (x,y,z),
        (x+width,y,z),
        (x+width,y+hight,z),
        (x,y+hight,z)
    ]

    

    glue = Glue(lcoords, depth)
    

    glue._finish_geometry()

    w1 = Waffle(lcoords = [(x,y-3*.1,z),(x+width,y-3*.1,z)],
                depth = depth,
                nLayers = 3,
                hightLayer = .1)

    w2 = Waffle(lcoords = [(x,y+hight,z),(x+width,y+hight,z)],
                depth = depth,
                nLayers = 3,
                hightLayer = .1)
    
    sandwich = Sandwich(w1, glue, w2)

    write_geo(sandwich.geometry, "SaaandwichClass")

def test17():
    s1, s2, s3, s4 = makeSomeSandwiches()

    lGeometries = [s1.geometry,s2.geometry,s3.geometry, s4.geometry]
    geometry =  mergeGeometries(lGeometries)
    write_geo(geometry, "UUhh")

def test18():
    s1, s2, s3, s4 = makeSomeSandwiches()
    rotorBlade = Blade(s1,s2,s3,s4)
    geometry = rotorBlade.geometry
    write_geo(geometry, "UUhh")

def test19():
    depth = 2
    width = 2

    hightGlue = .2
    hightLayer1 = 1
    hightLayer2 = 1

    totalHight = hightGlue+hightLayer1+hightLayer2
    nLayer1 = 2
    nLayer2 = nLayer1


    x = 0
    y = 0
    z = 0

    
    lcoordsGlue = [
        (x+.1, hightLayer1, z),
        (x+width-.1, hightLayer1, z),
        (x+width-.1, hightLayer1+hightGlue,z),
        (x+.1,hightLayer1+hightGlue,z)
    ]

    lcoordLayer1 = [
        (x,y,z),
        (x+width,y,z)
    ]

    lcoordLayer2 = [
        (x,y+hightLayer1+hightGlue,z),
        (x+width, y+hightLayer1+hightGlue,z)
    ]

    #=================================
    # FIRST SANDWICH
    #=================================

    glue = Glue(lcoordsGlue, depth)
    
    s = Sphere(.1, (0.4,1.15,-.5))
    glue.add_bubble(s)

    s = Sphere(.1, (0.6,0.95,-1.5))
    glue.add_bubble(s)

    glue._finish_geometry()

    w1 = Waffle(lcoords = lcoordLayer1,
                depth = depth,
                nLayers = nLayer1,
                hightLayer = hightLayer1/nLayer1)

    w2 = Waffle(lcoords = lcoordLayer2,
                depth = depth,
                nLayers = nLayer2,
                hightLayer = hightLayer2/nLayer2)
    
    sandwich1 = Sandwich(w1, glue, w2)
    write_geo(sandwich1.geometry, "Badabuff")


def test20():
    s1, s2, s3, s4 = makeSomeSandwiches2()
    rotorBlade = Blade(s1,s2,s3,s4)
    geometry = rotorBlade.geometry
    write_geo(geometry, "UUhh")