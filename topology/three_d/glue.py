from ..geo_utils.utils import write_geo
from .bubble import Sphere, clip

class Glue:
    # naming anti clockwise
    def __init__(self, lcoords, depth):
        """

        lines, loops etc everything anticlockwise, from outside perspektive


           c8__________c7
          / |         /|
         / c5        / |
        c4---------c3  c6
        | /         | /
        c1---------c2/
        
        lcoords = [c1, c2, c3, c4]

        faces naming:

        front: 1
        back: 2
        right:3
        left:4
        up:5
        down:6

        point = (x,y,z)

        y coordinate
         |
         |_____ x coordinate  
         /
        /
        z coordinate
        """
        self.height_up = lcoords[2][1]
        self.bubbles_in = []
        self.interior_lineLoops = []
        self.exterior_surfaces = []
        self.bubbles_on = {"up_interior":[],"down_interior":[]}    

        self.geo = {"free_point_id":9,"free_line_id":13,"free_lineLoop_id":7,"free_surface_id":1}

        # make points
        points = [{"id":i+1,"coords":lcoords[i],"lc":None} for i in range(4)]

        for i in range(4,8):
            points.append({"id":i+1,"coords":(lcoords[i%4][0],lcoords[i%4][1],lcoords[i%4][2]-depth), "lc":None})


        self.geo["points"]=points


        # make lines
        lines = []
        for i in range(4):
            lines.append({"id":i+1,"start":i+1,"end":(i+1)%4+1})
        
        for i in range(4,8):
            lines.append({"id":i+1, "start":i+1, "end":(i+1)%4+4+1})
        
        lines.append({"id":9,"start":3,"end":7})
        lines.append({"id":10,"start":2,"end":6})

        lines.append({"id":11,"start":4,"end":8})
        lines.append({"id":12,"start":1,"end":5})
            
        self.geo["lines"] = lines

        # make line loops
        line_loops = []

        line_loops.append({"id":1,"lines":[1,2,3,4]})
        line_loops.append({"id":2,"lines":[-5,-8,-7,-6]})
        line_loops.append({"id":3,"lines":[10,6,-9,-2]})
        line_loops.append({"id":4,"lines":[-12,-4,11,8]})
        line_loops.append({"id":5,"lines":[-3,9,7,-11]})
        line_loops.append({"id":6,"lines":[5,-10,-1,12]})

        self.geo["lineLoops"] = line_loops

        #self.faces = {i:{"lineLoop_main":i, "lineLoops_inside":[]} for i in range(1,7)}
        self.faces = {i: {"bubbles":[]} for i in range(1,7)}
        #self.faces = {"front":{},"back":{},"right":{},"left":{},"up":{},"down":{}}
        #self.faces["up"] = {"points":[4,3,7,8], "lines":[3,9,7,11], "lineLoops":[5]}
        #self.faces["down"] = {"points":[1,2,6,5],"lines":[1,10,5,12],"lineLoops":[6]}
        

    def _localize_bubble(self, bubble):
        return "inside", True 

    def add_bubble(self, bubble, loc, accuracy=.5):

        loc_, valid = self._localize_bubble(bubble)
        if valid:
            if loc == "inside":
                #self.bubbles_in.append(bubble)
                data = bubble.discretize(accuracy, self.geo["free_point_id"], self.geo["free_line_id"], self.geo["free_lineLoop_id"])
               
                self.geo["points"] += data["points"]
                self.geo["lines"] += data["lines"]
                lineLoops = []
                for triangle in data["lineLoops"]:

                    lineLoops.append({"id":triangle["id"], "lines":triangle["lines"]})

                self.geo["lineLoops"] += data["lineLoops"]
                self.interior_lineLoops += lineLoops

                self.geo["free_point_id"] += len(data["points"]) # TODO: correct this
                self.geo["free_line_id"] += len(data["lines"])
                self.geo["free_lineLoop_id"] += len(lineLoops)

            elif loc == "up":
                cut_height = self.height_up
                bubble_geo = bubble.discretize(accuracy, self.geo["free_point_id"], self.geo["free_line_id"], self.geo["free_lineLoop_id"])
                clip_bubble, intersectionLoop = clip(bubble_geo, cut_height)
                #print(clip_bubble["points"])
                self.geo["points"] += clip_bubble["points"]

                
                self.geo["lines"] += clip_bubble["lines"]
                self.geo["lineLoops"] += clip_bubble["lineLoops"]
                self.geo["lineLoops"].append(intersectionLoop)

                self.faces[5]["bubbles"].append({"clip_geo":clip_bubble,"intersectionLoop":intersectionLoop})
                
                self.geo["free_point_id"] = (max([p["id"] for p in clip_bubble["points"]])+1)
                self.geo["free_line_id"] = (max([l["id"] for l in clip_bubble["lines"]])+1)
                self.geo["free_lineLoop_id"]  = (max([ll["id"] for ll in clip_bubble["lineLoops"]])+2)
                
        else:
            print("Bubble could not be added")
            
    def _finish_geometry(self):
        self.geo["surfaces"] = []
        # Write exterior surfaces
        for face in self.faces:
            lintersectionLoops = []
            
            for bubble in self.faces[face]["bubbles"]:
                lintersectionLoops.append(bubble["intersectionLoop"]["id"])

            main_surface = {"id":face,"lineLoops":[face] + lintersectionLoops}
            self.geo["surfaces"].append(main_surface)
            self.exterior_surfaces.append(main_surface)

            for bubble in self.faces[face]["bubbles"]:
                llineLoops = bubble["clip_geo"]["lineLoops"]
                for lineLoop in llineLoops:
                    Id = lineLoop["id"]
                    surface = {"id":Id,"lineLoops":[Id]}
                    self.geo["surfaces"].append(surface)

                    self.exterior_surfaces.append(surface)
            

            

            #lineLoops = [self.faces[i]["lineLoop_main"]] + self.faces[i]["lineLoops_inside"]
            #surface = {"id":i, "lineLoops":lineLoops}
            #self.geo["surfaces"].append(surface)
            
        
        # Write interior surfaces
        for lineLoop in self.interior_lineLoops:
            self.geo["surfaces"].append({"id":lineLoop["id"], "lineLoops":[lineLoop["id"]]})
        
        # Write exterior surface loop
        self.geo["surfaceLoops"] = [{"id":1,"surfaces":[surface["id"] for surface in self.exterior_surfaces]}]
        #self.geo["surfaceLoops"] = [{"id":face, "surfaces":self.faces[face]["surfaces"]}for face in self.faces]
        # Write interior surface loop
        self.geo["surfaceLoops"].append({"id":2, "surfaces": [lineLoop["id"] for lineLoop in self.interior_lineLoops]})
        
        # Write volume
        self.geo["volumes"] = [{"id":1, "surfaceLoops":[1,2]}]

    def to_geo(self, file_name):
        self._finish_geometry()
        write_geo(self.geo, file_name)