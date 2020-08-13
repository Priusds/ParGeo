from abc import ABC, abstractmethod
from dolfin import *
import xml.etree.ElementTree as ET
import os
import numpy as np


class Bubble(ABC):
    @abstractmethod
    def trafo(self, theta, phi):
        """
            Centerd in (0,0,0)
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
        geo = discretize_reference(self.trafo, accuracy,point_id,line_id,lineLoop_id)
        new_points = []
        for point in geo["points"]:
            x,y,z = point["coords"]
            mx,my,mz = self.midpoint
            new_point = {"id":point["id"], "coords":[x+mx, y+my,z+mz], "lc":None}
            new_points.append(new_point)
        geo["points"] = new_points
        return geo

    def getBoundingBox(self):
        return [x + self.radius for x in self.midpoint] + [x - self.radius for x in self.midpoint]


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

def euclidean_projection(triangles_boundary,geo,remainingPoints, h, axis, section,P_):
    """
        changes remainingPoints
    """
    for p in P_:
        for pp in remainingPoints:
            if pp["id"] == p:
                coords = pp["coords"]
                pp["coords"][axis]

def minimization_projection(triangles_boundary,geo,remainingPoints, h, axis, section,P_):
    """
        geo = {"points":[...], "lines":[...], "lineLoops":[...]}
        triangles_boundary is a list containing indices of geo["lineLoops"] corresponding to triangles on the boundary
    """
    raise NotImplementedError

def new_clip(old_geo, clippingPlane, orientation, clippingRule):
    """
        clips old_geo at the clippingPlane, orientation tells us which halfspace is the one to be clipped away


    """
    # iterate over all triangles, check which one is in the good halfspace,
    #  bad halfspace or intersects with the clippingPlane
    for j,triangle in enumerate(geo["lineLoops"]):
        lline_id = triangle["lines"]
        lline_index = [old_geo["lines"][search_index("id",abs(l),geo["lines"])] for l in [l1,l2,l3]]
        points = list(set([l["start"] for l in lines] + [l["end"] for l in lines]))

def clip(geo, h, axis, section, project=euclidean_projection):
    """
    axis in [0,1,2]
    section > 0 or < 0
    if >0 cuts everything above h
    """
    if axis in [0,2]:
        raise NotImplementedError
    
    remainingPoints = []
    remainingLines = []
    remainingLineLoops = []
    clip_geo = {"points":remainingPoints, "lines":remainingLines, "lineLoops":remainingLineLoops}
    triangles_inside = []
    triangles_boundary = []

    L_in = []
    P_on = [] # list of point IDs of vertices of triangles which are on the boundary
    P_in = []
    P_ = [] # list of point IDs of vertices of triangles which are on the boundary and are outside
    lineLoop = [] # intersection lineLoop with cutting plane


    for j,triangle in enumerate(geo["lineLoops"]):
        
        l1, l2, l3 = triangle["lines"]
        lines = [geo["lines"][search_index("id",abs(l),geo["lines"])] for l in [l1,l2,l3]]
        points = list(set([l["start"] for l in lines] + [l["end"] for l in lines]))
        
        loc_sum,loc = locate(geo, points,h, axis, section) # 1 in, 0 mix, -1 out
        
        
        if loc_sum == 1: # IN
            triangles_inside.append(j)
            L_in += [abs(l) for l in triangle["lines"]]
            P_in += points

        elif loc_sum == 0: # IN/OUT
            # TODO: at this point decide between projecting the points or removing this triangle (project up)
            lines_on = [abs(l) for l in triangle["lines"]]
            triangles_boundary.append(j)
            P_on += points
            for i,inside in enumerate(loc):
                if not inside:
                    P_.append(points[i])
            
            for l in lines_on:
                p1 = geo["lines"][search_index("id",l,geo["lines"])]["start"]
                p2 = geo["lines"][search_index("id",l,geo["lines"])]["end"]

                if ((not locate_point(geo, p1, h)) and (not locate_point(geo, p2, h))):
                    lineLoop.append(l)
        else: # OUT
            pass
    for triangle in triangles_boundary:
        lines = geo["lineLoops"][triangle]["lines"] 
        lines = [abs(l) for l in lines]
        L_in += lines


    L_in = set(L_in)
    P_in = set(P_in)
    P_on = set(P_on)
    P_ = set(P_)

    L_tot = L_in


    remainingPoints += [geo["points"][search_index("id",i,geo["points"])] for i in list(P_in.union(P_on))]

    # Projection
    project(triangles_boundary,geo,remainingPoints,h,axis,section,P_) # changes remainingPoints
    
    remainingLines += [geo["lines"][search_index("id",i,geo["lines"])] for i in list(L_tot)]
    remainingLineLoops += [geo["lineLoops"][j] for j in (triangles_inside + triangles_boundary)]


    # Sort line Loop
    sorted_lineLoop = sort_lineLoop(geo, list(set(lineLoop)))
    free_lineLoop_id = max([i["id"] for i in geo["lineLoops"]])
    intersectionLoop = {"id":free_lineLoop_id+1,"lines":sorted_lineLoop}

    return clip_geo, intersectionLoop


def discretize_reference(trafo, accuracy, point_id=1,line_id=1,lineLoop_id=1):
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

    points = []
    for v in root[0][0]:
        data = v.attrib
        index = int(data["index"])
        coords = [data["x"], data["y"], data["z"]]
        points.append({"id":index + point_id,"coords":coords, "lc":None})
    surfaces = []
    lines = []
    line_index = 0
    for t in root[0][1]:
        data = t.attrib
        index = int(data["index"])
        surface = {"id":index + lineLoop_id, "lines":[]}
        data["v0"] = int(data["v0"])
        data["v1"] = int(data["v1"])
        data["v2"] = int(data["v2"])
        tlines = [[data["v0"] + point_id, data["v1"] + point_id], 
                [data["v1"]+ point_id, data["v2"]+ point_id], 
                [data["v2"]+ point_id, data["v0"]+ point_id]]
        
        for j in range(3):
            max_iter = len(lines)
            i = 0
            notFound = True
            while notFound and i < max_iter:
                line = [lines[i]["start"],lines[i]["end"]]
                if set(line) == set(tlines[j]):
                    notFound = False
                    if line == tlines[j]:
                        surface["lines"].append(lines[i]["id"])
                    else:
                        surface["lines"].append(-lines[i]["id"])
                    
                i += 1
            if notFound:
                line = {"id":line_index + line_id, "start":tlines[j][0], "end":tlines[j][1]}
                lines.append(line)
                surface["lines"].append(line["id"])
                line_index += 1
        assert len(surface["lines"]) == 3
        surfaces.append(surface)
    geo = {"points":points, "lines":lines, "lineLoops":surfaces}
    

    # TRANSFORM GEO
    
    #phis = []
    #thetas = []
    new_points = []
    for point in geo["points"]:
        x,y,z = point["coords"]
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
        new_point = {"id":point["id"], "coords":new_coords, "lc":None}
        new_points.append(new_point)

    geo["points"] = new_points
    return geo
    
def locate(geo, lpoint_id,h,axis,section):
    """
        lpoint_id is a list of three point Id, representing a triangle

        returns tuple = (t1, t2)
        t1 in {-1,0,1}
        t2 = loc
    """
    # TODO: point ID not -1
    lpoint_index = [search_index("id",Id,geo["points"]) for Id in lpoint_id]
    lpoint = [geo["points"][i] for i in lpoint_index]
    loc = [p["coords"][axis] < h for p in lpoint] if section > 0 else [p["coords"][axis] > h for p in lpoint]
    if sum(loc) == 3:
        # triangle is inside
        return 1,loc
    elif sum(loc) == 0:
        # triangle is outside
        return -1,loc
    else:
        #triangle is on border
        return 0,loc

def locate_point(geo, point, h):
    point_index = search_index("id", point,geo["points"])
    return geo["points"][point_index]["coords"][1] < h

def sort_lineLoop(geo, lines_id):
    lines = [geo["lines"][search_index("id",l,geo["lines"])] for l in lines_id]

    #lines = [geo["lines"][l-1] for l in lines_id]
    nlines = len(lines)
    free_lines = lines_id
    sorted_lines = [free_lines[0]]
    free_lines.pop(0)

    i = 0
    geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["start"] #geo["lines"][sorted_lines[-1]-1]["start"]

    while len(sorted_lines) < nlines and i < 1000:
        i+=1
        current = set([geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["start"],
                        geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["end"]])
        #current = set([geo["lines"][sorted_lines[-1]-1]["start"],geo["lines"][sorted_lines[-1]-1]["end"]])

        cond = True
        j = 0
        
        while cond and j < len(free_lines):
            free_line = set([geo["lines"][search_index("id", free_lines[j], geo["lines"])]["start"], geo["lines"][search_index("id", free_lines[j], geo["lines"])]["end"]])
            #free_line = set([geo["lines"][free_lines[j]-1]["start"], geo["lines"][free_lines[j]-1]["end"]])
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
        current_start = geo["lines"][search_index("id", abs(sorted_lines[i]), geo["lines"])]["start"]

        
        if i < nlines-1:
            next_start = geo["lines"][search_index("id", abs(sorted_lines[i+1]), geo["lines"])]["start"]
            next_end = geo["lines"][search_index("id", abs(sorted_lines[i+1]), geo["lines"])]["end"]
        else:
            next_start = geo["lines"][search_index("id", first_line_id, geo["lines"])]["start"]
            next_end = geo["lines"][search_index("id", first_line_id, geo["lines"])]["end"]

        
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