from .writeGeo import WriteGeo
from .star_geometry import  StarGeometry, Circle, Ellipse, Stellar, Rectangular, PolyRectangular
import os
import numpy as np
from dolfin import *
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as P
from matplotlib.collections import PatchCollection
import matplotlib
from shapely.geometry import Polygon, LinearRing



class Domain:
    # TODO: all constants should be relative
    DELTA = .018 # minimal distance between hole and Domain boundary
    EPSILON = .05 # minimal distance between hole and bounding box
    GAMMA = .02 # minimal distance between bounding boxes
    ETA = 2*EPSILON + GAMMA # minimal distance between holes
    NU = EPSILON - DELTA # minimal distance between boundary and bounding box
    ALLOWED_HOLES = {Circle, Stellar, Ellipse, Rectangular, PolyRectangular}

    # hole dependend mesh discretization
    RLAYERWIDTH = 0.02
    RALPHASHIFT = 0.5

    def __init__(self, boundary_hole):
        """
        boundary_hole is an obj. of class Hole
        """
        self.boundary = boundary_hole
        self.holes = []
        self.max_rads = [] # each hole with midpoint m is contained in a Ball(m,max_rad)
        self.bounding_boxes = []    # Just used in plot method in order to see the result of the cutout method.
        
    @classmethod
    def from_dict(cls, dict):
        return Domain(StarGeometry.from_dict(dict))

    def to_dict(self):
        dicts = {}
        dicts["boundary"] = self.boundary.to_dict()
        dicts["holes"] = {}
        for i, h in enumerate(self.holes):
            dicts["holes"][i] = h.to_dict()
        return dicts

    def add_hole(self, hole, do_test = True):
        """
        Addes hole to domain, if there are no collisions with the boundary and other holes, w.r.t. 
        DELTA, EPSILON, GAMMA, ETA.
        """
        if type(hole) not in Domain.ALLOWED_HOLES:
            print(type(hole))
            raise ValueError("Hole type not accepted.")

        if do_test and self._test(hole):
            self.holes.append(hole)
            return True
        elif do_test == False:
            self.holes.append(hole)
        else:
            return False


    def _test(self, hole, refs = 80):
        """
        Checks if the hole can be added to the Domain.
        
        Discretize all holes and check collisions beween polygons using shapely.
        """
        
        hole_ring = LinearRing(hole.discretize_hole(refs))

        discrete_boundary = self.boundary.discretize_hole(refs)
        max_rad = np.max(np.linalg.norm(np.array(discrete_boundary) - np.array(hole.midpoint), axis = 1))
        midpoint = np.array(hole.midpoint)


        # CHECK BOUNDARY COLLISION
        collision_with_boundary = True
        if self.boundary is Rectangular and self.boundary.angle == 0:
            
            r0, r1 = self.boundary.midpoint
            ax = self.boundary.width/2
            ay = self.boundary.height/2

            if r0 - ax + max_rad + Domain.DELTA <= midpoint[0] <= r0 + ax - max_rad - Domain.DELTA: 
                if r1 - ay + max_rad + Domain.DELTA <= midpoint[1] <= r1 + ay - max_rad - Domain.DELTA: 
                    collision_with_boundary = False
            print("First boundary test passed")
        else:
            boundary_poly = Polygon(discrete_boundary)
            boundary_ring = LinearRing(discrete_boundary)

            if not boundary_poly.contains(hole_ring):
                return False

            if hole_ring.distance(boundary_ring) <= Domain.DELTA:
                return False
            collision_with_boundary = False

        if collision_with_boundary:
            return False
        
        # CHECK HOLES COLLISION
        for rad, hole_ in zip(self.max_rads, self.holes):
            midpoint_ = hole_.midpoint
            distance = np.linalg.norm(np.array(midpoint_) - np.array(midpoint))
            if distance <= rad + max_rad + Domain.ETA:
                # Make accurate test
                h_ring = LinearRing(hole_.discretize_hole(refs))
                if h_ring.distance(hole_ring) <= Domain.ETA:
                    return False
        
        self.max_rads.append(max_rad)
        return True


    def _test_intersection_with_with_holes(self, hole, max_rad):
        midpoint = hole.midpoint
        


    def plot(self, refs_boundary = 60, refs_holes = 60, shownumbers = False):
        """
        Plots the actual domain.
        """
        holes_d = map(lambda x: np.array(x.discretize_hole(refs_holes)), self.holes + self.bounding_boxes)
        boundary_points = self.boundary.discretize_hole(refs_boundary)
        ax = plt.gca()
        eps = .2
        xlim_min = min(boundary_points, key = lambda x: x[0])[0] - eps
        xlim_max = max(boundary_points, key = lambda x: x[0])[0] + eps
        ylim_min = min(boundary_points, key = lambda x: x[1])[1] - eps
        ylim_max = max(boundary_points, key = lambda x: x[1])[1] + eps

        ax.set(xlim=(xlim_min, xlim_max), ylim=(ylim_min, ylim_max))

        patches = []
        for hole in holes_d:
            patches.append(P(hole))

        patches.append(P(boundary_points, linestyle = "-"))
        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

        ax.add_collection(p)

        if shownumbers:
            for s,h in enumerate(self.holes):
                mx,my = h.midpoint
                ax.text(mx,my,str(s))

        plt.autoscale(enable = True, tight = True, axis = 'both')


    def writeGeo(self, file_name, lc, refs_boundary, refs_holes,  agglomeration = False):
        file = WriteGeo(file_name)
        file.write_polygon(self.boundary.discretize_hole(refs_boundary), lc)
       

        if not agglomeration : 
            for hole in self.holes:
                d_hole = hole.discretize_hole(refs_holes)
                
                file.write_polygon(d_hole, lc)

            file.write_plane_surface(range(1, len(self.holes) + 2))
            file.close()

        else:
            for hole in self.holes:
                print(type(hole), "#########################################")
                d_hole_agg = hole.discretize_hole(refs_holes, rwidth = self.RLAYERWIDTH, ralphashift = self.RALPHASHIFT)
                file.write_polygon(d_hole_agg, lc)
            file.write_plane_surface(range(1, len(self.holes) + 2))

            index = len(self.holes) + 2
            for i, hole in enumerate(self.holes):
                d_hole = hole.discretize_hole(refs_holes)
                file.write_polygon(d_hole, lc) 
                file.write_plane_surface([i+2, index+i])   

            file.close()
        



    def safe_json(self, file_name):
        data = {}
        data["boundary"] = self.boundary.to_dict()
        for i, hole in enumerate(self.holes):
            data["hole" + str(i)] = hole.to_dict() #make sure no numpy array

        for i, box in enumerate(self.bounding_boxes):
            data["bounding_box" + str(i)] = box.to_dict() #make sure no numpy array

        with open(file_name + ".json", 'w') as fp:
            json.dump(data, fp, indent=2)


def load_json(file_name):
    with open(file_name + ".json", 'r') as fp:
        topo_json = json.load(fp)
   
    topo = Domain.from_dict(topo_json["boundary"])

    for key, hole in topo_json.items():
        if "bounding_box" in key:
            if not (hole["type"] == "rectangular" or hole["type"] == "poly_rectangular"):

                raise ValueError("Bounding box with type {0} is not rectangular".format(hole["type"]))
            new_rect = StarGeometry.from_dict(hole) 
            topo.bounding_boxes.append(new_rect)
        elif "hole" in key:
            topo.holes.append(StarGeometry.from_dict(hole))
    return topo

def dolfin_mesh(file_name):
    # creates .xml mesh from .geo mesh and returns dolfin mesh
    os.system("gmsh -2 " + file_name + ".geo")  # -v 0
    os.system("dolfin-convert " + file_name + ".msh " + file_name + ".xml")
    return Mesh(file_name + ".xml")