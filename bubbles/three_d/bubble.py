from abc import ABC, abstractmethod
import xml.etree.ElementTree as ET
import os
import numpy as np
from bubbles.three_d.utils import locate
from bubbles.three_d.clip import *


class Bubble(ABC):
    @abstractmethod
    def trafo(self, theta, phi):
        """
            Centered in (0,0,0)
            Returns the mapping from spheric coordinates to radius
            Implements the underlying geometric mapping from
          
            p = (theta, phi) in (0,pi) x (0, 2pi) --> F(p) in IR^3

        """
        pass

    @abstractmethod
    def discretize(self, args):
        """
            Returns an dictionary containging all infos needed for the .geo file
                - "points"  *   
                - "lines"   *
                - "lineLoops"   *
                - "surfaces"
            * these keys are obliged
        """
        pass

    @abstractmethod
    def getBoundingBox(self):
        """
            Returns [max_x, max_y, max_z, min_x, min_y, min_z]
        """
        pass

class Sphere(Bubble):
    def __init__(self, radius, midpoint):
        self.radius = radius
        self.midpoint = midpoint

    def trafo(self, theta, phi):
        return self.radius
    
    def discretize(self, accuracy, point_id=1,line_id=1,lineLoop_id=1):
        geometry = discretize_reference(self.trafo, accuracy,point_id,line_id,lineLoop_id)

        for point in geometry["points"]:
            x,y,z = geometry["points"][point]["coords"]
            mx,my,mz = self.midpoint
            new_coords = [x+mx, y+my,z+mz]
            geometry["points"][point]["coords"] = new_coords
     
        return geometry

    def getBoundingBox(self):
        return [x + self.radius for x in self.midpoint] + [x - self.radius for x in self.midpoint]

class Ellipse(Bubble):
    def __init__(self, radii, midpoint):
        self.radii = radii
        self.midpoint = midpoint

    def trafo(self, theta, phi):
        """
            theta in (0, pi)
            phi in (0, 2pi)
        """
        theta = theta - np.pi/2
        a,b,c = self.radii
        return (a*b*c)/np.sqrt(c**2*(b**2*np.cos(phi)**2 + a**2*np.sin(phi)**2)*np.cos(theta)**2+a**2*b**2*np.sin(theta)**2)
    
    def discretize(self, accuracy, point_id=1,line_id=1,lineLoop_id=1):
        geometry = discretize_reference(self.trafo, accuracy,point_id,line_id,lineLoop_id)

        for point in geometry["points"]:
            x,y,z = geometry["points"][point]["coords"]
            mx,my,mz = self.midpoint
            new_coords = [x+mx, y+my,z+mz]
            geometry["points"][point]["coords"] = new_coords
     
        return geometry

    def getBoundingBox(self):
        return [x + max(self.radii) for x in self.midpoint] + [x - max(self.radii) for x in self.midpoint]

class StarBubble(Bubble):
    def __init__(self, size, midpoint, coefficients):
        self.midpoint = midpoint
        self.size = size
        self.coefficients = coefficients

    def trafo(self, theta, phi):
        """
            @param p = (theta, phi)   theta in (0,pi)
                                      phi   in (0,2pi)


                  p[0] is an array of thetas etc

                  p in IR^(2 x N)

        """
        p = (theta,phi)
        res = np.zeros(3)
        a,b,c = 1,1,1
        res[0] = (1+0.2*np.sin(3*p[0]) + .3*np.cos(5*p[1]))*a * np.sin(p[0]) * np.cos(p[1]) # (1+0.2*np.sin(3*p[0]) + .3*np.cos(5*p[1]))*
        res[1] = b*(1+0.2*np.sin(3*p[0]) + .3*np.cos(5*p[1])) * np.sin(p[0]) * np.sin(p[1])
        res[2] = c * np.cos(p[0])

        return np.linalg.norm(res)*.1
    
    def discretize(self, accuracy=.1, point_id=1,line_id=1,lineLoop_id=1):
        geo = discretize_reference(self.trafo, accuracy,point_id,line_id,lineLoop_id)
        new_points = []
        for point in geo["points"]:
            x,y,z = point["coords"]
            mx,my,mz = self.midpoint
            new_point = {"id":point["id"], "coords":[x+mx, y+my,z+mz], "lc":None}
            new_points.append(new_point)
        geo["points"] = new_points
        return geo

    def getBoundingBox(self, nsample=40):
        # sample based
        mx,my,mz = self.midpoint
        ltheta = np.linspace(0,np.pi,nsample)
        lphi = np.linspace(0,2*np.pi,nsample)
        points = [_targetPoint(theta,phi,self.trafo(theta, phi))  for theta, phi in zip(ltheta,lphi)]
        points = [(p[0]+mx, p[1]+my, p[2]+mz) for p in points]

        max_x = max(points, key=lambda x: x[0])[0]
        max_y = max(points, key=lambda x: x[1])[1]
        max_z = max(points, key=lambda x: x[2])[2]

        min_x = min(points, key=lambda x: x[0])[0]
        min_y = min(points, key=lambda x: x[1])[1]
        min_z = min(points, key=lambda x: x[2])[2]

        return [max_x, max_y, max_z, min_x, min_y, min_z]

def _targetPoint(theta, phi, radius):
    x = radius*np.sin(theta)*np.cos(phi)
    y = radius*np.sin(theta)*np.sin(phi)
    z = radius*np.cos(theta)
    return (x,y,z)


def clip(geometry, clippingPlane, pointProjection,clippingRule, args = None):
    """
        clips old_geo at the clippingPlane, orientation tells us which halfspace is the one to be clipped away


    """

    # iterate over all triangles, check which one is in the good halfspace,
    #  bad halfspace or intersects with the clippingPlane
    if args is not None : 
        clipper = clippingRule(pointProjection(clippingPlane, args))
    else:
        clipper = clippingRule(pointProjection(clippingPlane))

    normalVector, supportVector = clippingPlane


    ltriangle_inside = set()
    lline_inside = set()
    lpoint_inside = set()

    ltriangle_boundary = set()
    lline_boundary = set()
    lpoint_boundary = set()

    lpoint_outside = set()
    lineLoop = set()

    for triangle, lines in geometry["lineLoops"].items():
        lines_id = [abs(l) for l in lines]
        points = {geometry["lines"][l][0] for l in lines_id}.union({geometry["lines"][l][1] for l in lines_id})
        points_locations = {point:locate(geometry["points"][point]["coords"], normalVector, supportVector)  for point in points}
        s = sum([_ for _ in points_locations.values()])
        #print(points)
        if s == 3:
            # Triangle is inside
            ltriangle_inside.add(triangle)
            lline_inside = lline_inside.union(lines_id)
            lpoint_inside = lpoint_inside.union(points)

        elif 1 <= s <= 2:
            # Triangle is on the boundary
            ltriangle_boundary.add(triangle)
            lline_boundary = lline_boundary.union(lines_id)
            lpoint_boundary = lpoint_boundary.union(points)

            for l_id in lines_id:
                start,end = geometry["lines"][l_id]
                if (not points_locations[start] and  not points_locations[end]):
                    lineLoop.add(l_id)


        else:
            # Triangle is outside
            pass

    ltriangle_boundary, lline_boundary, new_points = clipper.apply(ltriangle_boundary, lline_boundary, lpoint_boundary, geometry) 

    lpoint_inside = lpoint_inside.difference(lpoint_boundary)
    lline_inside = lline_inside.difference(lline_boundary)
    clip_geometry = dict()
    clip_geometry["points"] = {**{Id:geometry["points"][Id] for Id in lpoint_inside},**new_points}
    clip_geometry["lines"] = {**{Id:geometry["lines"][Id] for Id in lline_inside}, **{Id:geometry["lines"][Id] for Id in lline_boundary}}
    clip_geometry["lineLoops"] = {**{Id:geometry["lineLoops"][Id] for Id in ltriangle_inside}, **{Id:geometry["lineLoops"][Id] for Id in ltriangle_boundary}}
    
    assert len(lineLoop) >= 3
    
    sorted_lineLoop = sort_lineLoop(geometry, list(lineLoop))
    free_lineLoop_id = max([Id+1 for Id in geometry["lineLoops"]])
    intersectionLoop = (free_lineLoop_id,sorted_lineLoop) 

    return clip_geometry, intersectionLoop  
    

def get_key(val,my_dict): 
    for key, value in my_dict.items(): 
        if val == value: 
            return key
        elif val == (value[1],value[0]):
            return key
    return None
def discretize_reference(trafo, accuracy, point_id=1,line_id=1,lineLoop_id=1):
    # TODO: correct geo to geometry
    # 1: Write geo
    # 2: Write xml
    # 3: Load xml TODO: speed up
    # 4: transform
    # return geomerty: dict containing points,lines,lineLoops

    # 1: Write geo
    reference_sphere = open("sphere.geo", "r")

    finer_sphere = open("sphere_finer.geo", "w")

    txt = []
    lc = accuracy
    for i,line in enumerate(reference_sphere):

        if i == 0:

            new_line = "h = " + str(lc) + ";\n"
            finer_sphere.write(new_line)
        else:
            finer_sphere.write(line)
    reference_sphere.close()
    finer_sphere.close()
    #  2:   WRITE XML

    os.system("gmsh -2 " + "sphere_finer" + ".geo")  # -v 0
    os.system("dolfin-convert " + "sphere_finer" + ".msh " + "sphere_finer" + ".xml")
    
    # 3:    LOAD XML 
    tree = ET.parse('sphere_finer.xml')

    root = tree.getroot()
    dpoints = {int(v.attrib["index"])+point_id : {"coords":[v.attrib["x"], v.attrib["y"], v.attrib["z"]],"lc":None} for v in root[0][0]}
    surfaces = []
    lineLoops = dict()
    dlines = dict()
    line_index = 0
    for t in root[0][1]:
        data = t.attrib
        index = int(data["index"])
        surface = {"id":index + lineLoop_id, "lines":[]}
        data["v0"] = int(data["v0"])
        data["v1"] = int(data["v1"])
        data["v2"] = int(data["v2"])
        tlines = [(data["v0"] + point_id, data["v1"] + point_id), 
                (data["v1"]+ point_id, data["v2"]+ point_id), 
                (data["v2"]+ point_id, data["v0"]+ point_id)]
        lineLooplines = []
        
        for line in tlines:
            Id = get_key(line, dlines)
            if Id is None:
                dlines[line_index+line_id] = line
                Id = line_index+line_id
                line_index += 1
            else:
                if dlines[Id][0] != line[0]:
                    Id = -Id
            lineLooplines.append(Id)
        lineLoops[index+lineLoop_id] = lineLooplines
        assert len(lineLooplines) == 3
            
    geo = {"points":dpoints, "lines":dlines, "lineLoops":lineLoops}
    
    # TRANSFORM GEO
    
    for point in geo["points"]:
        x,y,z = geo["points"][point]["coords"]
        x = float(x)
        y = float(y)
        z = float(z)
        phi = np.arctan2(y,x) + 2*np.pi
        #phis.append(phi)
        theta = np.arccos(z)
        #thetas.append(theta)
        
        new_radius = trafo(phi, theta)
        new_x = new_radius*x
        new_y = new_radius*y
        new_z = new_radius*z
        new_coords = [new_x, new_y, new_z]
        geo["points"][point]["coords"] = new_coords
        
    return geo



def sort_lineLoop(geometry, lines_id):
    nlines = len(lines_id)
    free_lines = lines_id.copy()

    sorted_lines = [free_lines[0]]
    free_lines.pop(0)

    maxIter = 1000
    i=0
    while len(free_lines)>0 and i<maxIter:
        head_ID = sorted_lines[0]
        tail_ID = sorted_lines[-1]

        head = set(geometry["lines"][head_ID])
        tail = set(geometry["lines"][tail_ID])


        cond = True
        j = 0
        while cond and j < len(free_lines):
            candidate = set(geometry["lines"][free_lines[j]])

            if len(candidate.intersection(tail)) > 0:
                sorted_lines.append(free_lines[j])
                free_lines.remove(free_lines[j])
                cond = False
            elif len(candidate.intersection(head)) > 0:
                sorted_lines = [free_lines[j]] + sorted_lines
                free_lines.remove(free_lines[j])
                cond = False
            j += 1
        
        i+=1

    i = 0
    while len(sorted_lines) < nlines and i < 0:
        i+=1
        current = set([geometry["lines"][sorted_lines[-1]][0],geometry["lines"][sorted_lines[-1]][1]])
        cond = True
        j = 0
        
        while cond and j < len(free_lines):
            free_line = set([geometry["lines"][free_lines[j]][0], geometry["lines"][free_lines[j]][1]])
            if len(free_line.intersection(current)) > 0:
                sorted_lines.append(free_lines[j])
                free_lines.remove(free_lines[j])
                cond = False
            j += 1
    if len(sorted_lines) != nlines:
        raise ValueError("Problem with clipping: no sorting for lineLoop")

    # Correct line sign for lineLoop
    first_line_id = abs(sorted_lines[0])
    for i in range(nlines):
        current_start = geometry["lines"][abs(sorted_lines[i])][0]
        
        if i < nlines-1:
            next_start = geometry["lines"][abs(sorted_lines[i+1])][0]
            next_end = geometry["lines"][abs(sorted_lines[i+1])][1]
        else:
            next_start = geometry["lines"][first_line_id][0]
            next_end = geometry["lines"][first_line_id][1]

        
        if current_start in {next_start,next_end}:
            sorted_lines[i] = -sorted_lines[i]
        
    return sorted_lines

def search_index(key, value, list_of_dicts):
    l = len(list_of_dicts)
    i = 0
    while i < l:
        if list_of_dicts[i][key] == value:
            return i
        i += 1
    raise ValueError("Did not find ", value, " of ", key)