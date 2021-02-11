from topology.two_d.topology import Topology
from topology.two_d.hole import Stellar, Circle
from topology.geo_utils.utils import write_geo

import random as rd
from dolfin import *
import os
import matplotlib.pyplot as plt
import numpy as np
import json

def dolfin_mesh(file_name):
    # creates .xml mesh from .geo mesh and returns dolfin mesh
    os.system("gmsh -2 " + file_name + ".geo")  # -v 0
    os.system("dolfin-convert " + file_name + ".msh " + file_name + ".xml")
    return Mesh(file_name + ".xml")



def couldNotProcessBug(fname):
    width = 1.
    height = 1.
    domain = {"x0" : 0, 
              "y0" : 0, 
              "w"  : width,
              "h"  : height
    }
    x0, y0 = (domain["x0"], domain["y0"]) 

    N = 2
    nHoles = 100
    filled = True

    rect_out = [(x0,y0),(x0 + width, y0 + height)]

    rects_in = []   

    # local placement of rectangles   s1 free  s2 rect  s3 free
    s1, s2, s3 = 1./4, 1./2, 1./4    # must sum to 1
        
    w_loc = width / N   # local width  of a N equidistant subdivided rectangles from Domain one
    h_loc = height / N  # local height of a N equidistant subdivided rectangles from Domain one

    xshift = s1*w_loc  # midpoint offset for rectangle in x direction
    yshift = s1*h_loc  # midpoint offset for rectangle in y direction

    w = s2 * w_loc  # width of subrectangle
    h = s2 * h_loc  # height of subrectangle

    for i in range( N ) : 
        for j in range( N ) : 
            # midpoint of new sub rectangle
            M = xshift + i * w_loc,  yshift + j * h_loc
            rect = [(M[0],M[1]),(M[0]+w,M[1]+h)]
            rects_in.append(rect)           
        
    discrete_Holes = []

    refs = 40      # refs for the interior inclusions
    periodic_boundary = False



    M = 100 


    for m in range(M):

        print(" =================== START NEW TOPOLOGY ===================")

        topo = Topology(rect_out, rects_in, periodic_boundary= periodic_boundary)


        hNa = 0    # attempted hole number
        ddict = {} # debug storage dict


        nCounter = 0
        #for k in range(100):
        while nCounter < nHoles:
                    
            midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
            radius = 0.015 * min(width, height)  #rd.uniform(0,.1)
            hole = Circle(midpoint_s,radius) 

            ddict[hNa] = [hole.to_dict() , refs, False]

            hNa +=1


            with open(fname + ".json", 'w') as fp:
                json.dump(ddict, fp, indent=2)

                        
            if topo.add_hole(hole,refs=refs): 
                nCounter += 1

                ddict[hNa-1][2] = True  # hole has been added

def CouldNotProcessBugMinimal(verbose = False):
    width = 1.
    height = 1.
    domain = {"x0" : 0, 
              "y0" : 0, 
              "w"  : width,
              "h"  : height
    }
    x0, y0 = (domain["x0"], domain["y0"]) 

    N = 2
    filled = True

    rect_out = [(x0,y0),(x0 + width, y0 + height)]

    rects_in = []   

    # local placement of rectangles   s1 free  s2 rect  s3 free
    s1, s2, s3 = 1./4, 1./2, 1./4    # must sum to 1
        
    w_loc = width / N   # local width  of a N equidistant subdivided rectangles from Domain one
    h_loc = height / N  # local height of a N equidistant subdivided rectangles from Domain one

    xshift = s1*w_loc  # midpoint offset for rectangle in x direction
    yshift = s1*h_loc  # midpoint offset for rectangle in y direction

    w = s2 * w_loc  # width of subrectangle
    h = s2 * h_loc  # height of subrectangle

    for i in range( N ) : 
        for j in range( N ) : 
            # midpoint of new sub rectangle
            M = xshift + i * w_loc,  yshift + j * h_loc
            rect = [(M[0],M[1]),(M[0]+w,M[1]+h)]
            rects_in.append(rect)           
  

    refsNew = 40      # the number of refs does not change the appearance of the bug
    periodic_boundary = False

    topo = Topology(rect_out, rects_in, periodic_boundary= periodic_boundary)


    if verbose:
        plt.plot( [domain["x0"], domain["x0"]+domain["w"], domain["x0"]+domain["w"], domain["x0"], domain["x0"]],
                    [domain["y0"], domain["y0"], domain["y0"]+domain["h"], domain["y0"]+domain["h"], domain["y0"]],
                    'k--'
            )

        for rect in rects_in : 
            A, B = rect[0], rect[1]

            plt.plot([A[0], B[0], B[0], A[0], A[0]], [A[1],A[1],B[1],B[1],A[1]], 'k--', linewidth = 0.5)


    with open("bugHole" + ".json", 'r') as fp:
        holes_json = json.load(fp)
        k = '0'
        hole_data, refs = holes_json[k]

        M , r = hole_data['midpoint'], hole_data['radius']

        if verbose : 
            s = 2*pi / refs
            x = [M[0] + r * np.cos( n * s * 2*np.pi) for n in range(refsNew)]
            y = [M[1] + r * np.sin( n * s * 2*np.pi) for n in range(refsNew)]
            plt.plot(x,y, 'r.', ms = 1)
            plt.show()




        hole = Circle(M, r)
        topo.add_hole(hole,refs=refsNew)

      

                













def ExtractProcessBug():
    width = 1.
    height = 1.
    domain = {"x0" : 0, 
              "y0" : 0, 
              "w"  : width,
              "h"  : height
    }
    x0, y0 = (domain["x0"], domain["y0"]) 

    N = 2
    nHoles = 100
    filled = True

    rect_out = [(x0,y0),(x0 + width, y0 + height)]

    rects_in = []   

    # local placement of rectangles   s1 free  s2 rect  s3 free
    s1, s2, s3 = 1./4, 1./2, 1./4    # must sum to 1
        
    w_loc = width / N   # local width  of a N equidistant subdivided rectangles from Domain one
    h_loc = height / N  # local height of a N equidistant subdivided rectangles from Domain one

    xshift = s1*w_loc  # midpoint offset for rectangle in x direction
    yshift = s1*h_loc  # midpoint offset for rectangle in y direction

    w = s2 * w_loc  # width of subrectangle
    h = s2 * h_loc  # height of subrectangle

    for i in range( N ) : 
        for j in range( N ) : 
            # midpoint of new sub rectangle
            M = xshift + i * w_loc,  yshift + j * h_loc
            rect = [(M[0],M[1]),(M[0]+w,M[1]+h)]
            rects_in.append(rect)           
        
    discrete_Holes = []

    refs = 40      # refs for the interior inclusions
    periodic_boundary = False


    topo = Topology(rect_out, rects_in, periodic_boundary= periodic_boundary)

    file_name = "notProcessBug"

    with open(file_name + ".json", 'r') as fp:
        holes_json = json.load(fp)

        #print(holes_json)

        k = 0

        while str(k) in holes_json:
            hole_data, refs, added = holes_json[str(k)]

            print(str(k), " : ", holes_json[str(k)], "\n")

            hole = Circle(hole_data['midpoint'], hole_data['radius'])
            if str(k) == '36':
                bugHole = {}
                bugHole[0] = [hole.to_dict() , refs]

                with open("bugHole" + ".json", 'w') as fp:
                    json.dump(bugHole, fp, indent=2)


                topo.add_hole(hole,refs=refs)

            k+=1





def debugExamples(fname, loadbug = False):
    width = 1.
    height = 1.
    domain = {"x0" : 0, 
              "y0" : 0, 
              "w"  : width,
              "h"  : height
    }
    x0, y0 = (domain["x0"], domain["y0"]) 

    N = 2
    nHoles = 10
    filled = True

    def createTopology(x0, y0, periodic_boundary=False):

        rect_out = [(x0,y0),(x0 + width, y0 + height)]

        rects_in = []   

        # local placement of rectangles   s1 free  s2 rect  s3 free
        s1, s2, s3 = 1./4, 1./2, 1./4    # must sum to 1
        
        w_loc = width / N   # local width  of a N equidistant subdivided rectangles from Domain one
        h_loc = height / N  # local height of a N equidistant subdivided rectangles from Domain one

        xshift = s1*w_loc  # midpoint offset for rectangle in x direction
        yshift = s1*h_loc  # midpoint offset for rectangle in y direction

        w = s2 * w_loc  # width of subrectangle
        h = s2 * h_loc  # height of subrectangle

        for i in range( N ) : 
            for j in range( N ) : 
                # midpoint of new sub rectangle
                M = xshift + i * w_loc,  yshift + j * h_loc
                rect = [(M[0],M[1]),(M[0]+w,M[1]+h)]
                rects_in.append(rect)           
        


        discrete_Holes = []

    
        refs = 40      # refs for the interior inclusions


        if not loadbug : 
            # Set periodic_boundary, True or False for outer rectangle
            topo = Topology(rect_out, rects_in, periodic_boundary= periodic_boundary)

            nCounter = 0
            #for k in range(100):
            while nCounter < nHoles:
                    
                midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
                radius = 0.015 * min(width, height)  #rd.uniform(0,.1)
                stellar_hole = Circle(midpoint_s,radius) 
                        
                if topo.add_hole(stellar_hole,refs=refs): 
                    discrete_Holes.append(stellar_hole.discretize_hole(refs))
                    nCounter += 1

                    # make a random composite material via filled setted to true option
            topo.safe_json(fname)
                    
            # write geo file
            topo.write_geo(fname, filled = filled, lc_in=.025, lc_out=.025)
            print("Geo file written : ", fname)

        else : 
            topo = Topology.load_topology("bug")
        

        
    createTopology(x0, y0, periodic_boundary = False)

    mesh = dolfin_mesh(fname)
    #plot(mesh)
    #plt.show()



def main():
    CouldNotProcessBugMinimal(True)
    #ExtractProcessBug()
    #couldNotProcessBug("notProcessBug")

    #fname = "sample_20_Nsub2_Nholes100"
    #dolfin_mesh(fname)                  # this fails since the geo file seems to be corrupt
    #for k in range(10):
    #    debugExamples("bug", loadbug = True)



if __name__ == "__main__":
        main()  # sys.argv)