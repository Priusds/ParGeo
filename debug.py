from topology.two_d.topology import Topology
from topology.two_d.hole import Stellar, Circle
from topology.geo_utils.utils import write_geo

import random as rd
from dolfin import *
import os
import matplotlib.pyplot as plt

def dolfin_mesh(file_name):
    # creates .xml mesh from .geo mesh and returns dolfin mesh
    os.system("gmsh -2 " + file_name + ".geo")  # -v 0
    os.system("dolfin-convert " + file_name + ".msh " + file_name + ".xml")
    return Mesh(file_name + ".xml")

def debugExamples(fname):
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

        return discrete_Holes,  rects_in
    
    discrete_Holes, rects_in = createTopology(x0, y0, periodic_boundary = False)

    mesh = dolfin_mesh(fname)
    plot(mesh)
    plt.show()



def main():

    #fname = "sample_20_Nsub2_Nholes100"
    #dolfin_mesh(fname)                  # this fails since the geo file seems to be corrupt

    debugExamples("test")



if __name__ == "__main__":
        main()  # sys.argv)