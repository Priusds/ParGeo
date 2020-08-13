from ..geo_utils.utils import write_geo
from .bubble import Sphere, clip, StarBubble

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
        self.max_x = lcoords[1][0]
        self.max_y = lcoords[2][1]
        self.max_z = 0
        self.min_x = lcoords[0][0]
        self.min_y = lcoords[0][1]
        self.min_z = -depth
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
        
        self.geometry = {"points":{p["id"]: {"coords":p["coords"], "lc":None} for p in points}}
        print(self.geo)
        print(self.geometry)

    def _localize_bubble(self, bubble):
        # Assume midpoint is inside of the box and the bubble is small enough not to intersect with opposite facets
        max_x,max_y,max_z,min_x,min_y,min_z = bubble.getBoundingBox()
        valid = True
        loc = ""
        if max_x >= self.max_x:
            #intersection with right
            if max_y >= self.max_y:
                #intersection with up
                if max_z >= self.max_z:
                    #intersection with front
                    # cornercase
                    loc += "right_up_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    # cornercase
                    loc += "right_up_back"
                    valid = False
                else:
                    loc += "right_up"
                    valid = False
                #no intersection with right/left
            elif min_y <= self.min_y:
                #intersection with down
                if max_z >= self.max_z:
                    #intersection with front
                    # cornercase
                    loc += "right_down_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    # cornercase
                    loc += "right_down_back"
                    valid = False
                else:
                    #no intersection with front/back
                    loc += "right_down"
                    valid = False
            else:
                #no intersection with up/down
                if max_z >= self.max_z:
                    #intersection with front
                    loc += "right_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    loc += "right_back"
                    valid = False
                else:
                    #no intersection with front/back
                    loc += "right"
                    valid = True
                
        elif min_x <= self.min_x:
            #intersection with left
            if max_y >= self.max_y:
                #intersection with up
                if max_z >= self.max_z:
                    #intersection with front
                    #cornercase
                    loc += "left_up_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    #cornercase
                    loc += "left_up_back"
                    valid = False
                else:
                    #no intersection with front/back
                    valid = False

                #intersection with up
            elif min_y <= self.min_y:
                #intersection with down
                if max_z >= self.max_z:
                    #intersection with front
                    #cornercase
                    loc += "left_down_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    #cornercase
                    loc += "left_down_back"
                    valid = False
                else:
                    #no intersection with front/back
                    valid = False

            else:
                #no intersection with up/down
                if max_z >= self.max_z:
                    #intersection with front
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    valid = False
                else:
                    #no intersection with front/back
                    valid = True
                    loc += "left"
        else:
            # no intersection with left/right
            if max_y >= self.max_y:
                #intersection with up
                if max_z >= self.max_z:
                    #intersection with front
                    loc += "up_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    loc += "up_back"
                    valid = False
                else:
                    #no intersection with front/back
                    loc = "up"
            elif min_y <= self.min_y:
                #intersection with down
                if max_z >= self.max_z:
                    #intersection with front
                    loc += "down_front"
                    valid = False
                elif min_z <= self.min_z:
                    #intersection with back
                    loc += "down_back"
                    valid = False
                else:
                    #no intersection with front/back
                    loc = "down"
            else:
                #no intersection with up/down
                if max_z >= self.max_z:
                    #intersection with front
                    loc += "front"
                elif min_z <= self.min_z:
                    #intersection with back
                    loc += "back"
                else:
                    #no intersection with anyside
                    loc += "inside"
        
        return loc, valid

    def add_bubble(self, bubble, accuracy=.5):

        loc, valid = self._localize_bubble(bubble)
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
                clip_bubble, intersectionLoop = clip(bubble_geo, cut_height, 1,1)
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
                raise NotImplementedError(loc)
        else:
            print("Bubble could not be added", loc, valid)
            
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