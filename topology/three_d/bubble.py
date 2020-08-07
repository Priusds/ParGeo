from abc import ABC, abstractmethod
from dolfin import *
import xml.etree.ElementTree as ET
import os
import numpy as np


class Bubble(ABC):
    @abstractmethod
    def trafo(self, theta, phi):
        """
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

class StarBubble(Bubble):
    def __init__(self, midpoint, size, coefficients):
        self.midpoint = midpoint
        self.size = size
        self.coefficients = coefficients

    def trafo(self, theta, phi):
        raise NotImplementedError
    
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

def clip(geo, h):
    finalPoints = []
    finalLines = []
    finalLineLoops = []
    clip_geo = {"points":finalPoints, "lines":finalLines, "lineLoops":finalLineLoops}
    triangles_total = geo["lineLoops"]
    triangles_inside = []
    triangles_outside = []
    triangles_inout = []
    L_out = []
    L_on = []
    L_in = []
    P_out = []
    P_on = []
    P_in = []
    P_ = []
    lineLoop = []
    i = 1

    for triangle in geo["lineLoops"]:
        
        l1, l2, l3 = triangle["lines"]
        lines = [geo["lines"][search_index("id",abs(l),geo["lines"])] for l in [l1,l2,l3]]
        #print([geo["lines"][l-1] for l in [l1,l2,l3]])
        #points = list(set([geo["lines"][abs(l)-1]["start"] for l in triangle["lines"]] + [geo["lines"][abs(l)-1]["end"] for l in triangle["lines"]]))
        points = list(set([l["start"] for l in lines] + [l["end"] for l in lines]))
        
        loc_sum,loc = locate(geo, points, h) # 1 in, 0 mix, -1 out
        
        if loc_sum == 1: # IN
            triangles_inside.append(triangle["id"])
            L_in += [abs(l) for l in triangle["lines"]]
            P_in += points
        elif loc_sum == -1: # OUT
            triangles_outside.append(triangle["id"])
            P_out += points
            L_out += [abs(l) for l in triangle["lines"]]
        else: # IN/OUT
            lines_on = [abs(l) for l in triangle["lines"]]
            L_on += lines_on
            triangles_inout.append(triangle["id"])
            P_on += points
            for i,boolean in enumerate(loc):
                if not boolean:
                    P_.append(points[i])
            
            for l in lines_on:
                p1 = geo["lines"][search_index("id",l,geo["lines"])]["start"]
                p2 = geo["lines"][search_index("id",l,geo["lines"])]["end"]
                #p1 = geo["lines"][l-1]["start"]
                #p2 = geo["lines"][l-1]["end"]
                if ((not locate_point(geo, p1, h)) and (not locate_point(geo, p2, h))):
                    lineLoop.append(l)


    for triangle in triangles_inout:
        triangle_index = search_index("id",triangle,geo["lineLoops"])
        lines = geo["lineLoops"][triangle_index]["lines"] # lines = geo["lineLoops"][triangle -1 ]["lines"] 
        lines = [abs(l) for l in lines]
        L_in += lines

    L_out = set(L_out)
    L_on = set(L_on)
    L_in = set(L_in)
    P_out = set(P_out)
    P_in = set(P_in)
    P_on = set(P_on)

    L_tot = L_in # .union(L_out.intersection(L_on))


    finalPoints += [geo["points"][search_index("id",i,geo["points"])] for i in list(P_in.union(P_on))]
    #finalPoints += [geo["points"][i-1] for i in list(P_in.union(P_on))]

    # Process finalPoints # TODO: Write application
    for p in P_:
        for pp in finalPoints:
            if pp["id"] == p:
                coords = pp["coords"]
                pp["coords"] = [coords[0], h ,coords[2]]
    
    finalLines += [geo["lines"][search_index("id",i,geo["lines"])] for i in list(L_tot)]
    #finalLines += [geo["lines"][i-1] for i in list(L_tot)]
    finalLineLoops += [geo["lineLoops"][search_index("id",i,geo["lineLoops"])] for i in (triangles_inside + triangles_inout)]
    #finalLineLoops += [geo["lineLoops"][i-1] for i in (triangles_inside + triangles_inout)]

    # Sort line Loop
    sorted_lineLoop = sort_lineLoop(geo, list(set(lineLoop)))
    free_lineLoop_id = max([i["id"] for i in geo["lineLoops"]])
    intersectionLoop = {"id":free_lineLoop_id+1,"lines":sorted_lineLoop}
    #clip_geo["lineLoops"].append(intersectionLoop)
    #clip_geo["surfaces"] = [{"id":ll["id"], "lineLoops":ll["id"]} for ll in clip_geo["lineLoops"]]
    #write_geo(clip_geo, "test_bubbly")
    return clip_geo, intersectionLoop


def discretize_reference(trafo, accuracy, point_id=1,line_id=1,lineLoop_id=1):
    # 1: Write geo
    # 2: Write xml
    # 3: Load xml
    # 4: transform
    # return geomerty: dict containing points,lines,lineLoops

    # 1:
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
    #print(min(phis), max(phis))
    #print(min(thetas), max(thetas))
    geo["points"] = new_points
    return geo
    
def locate(geo, lpoint_id,h=-.06):
    # TODO: point ID not -1
    #self.geo["point_id"], self.geo["line_id"], self.geo["lineLoop_id"]print(geo["points"])
    lpoint_index = [search_index("id",Id,geo["points"]) for Id in lpoint_id]
    lpoint = [geo["points"][i] for i in lpoint_index]
    loc = [p["coords"][1] < h for p in lpoint] # false if out
    if sum(loc) == 3:
        return 1,loc
    elif sum(loc) == 0:
        return -1,loc
    else:
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