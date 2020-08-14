from ..geo_utils.utils import write_geo
from .bubble import Sphere, clip, StarBubble
from .clip import EasyClip, EuclideanProjection
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

        self.freePointId = 9
        self.freeLineId = 13
        self.freeLineLoopId = 7
        self.freeSurfaceId = 1

        # make points
        points = [{"id":i+1,"coords":lcoords[i],"lc":None} for i in range(4)]

        for i in range(4,8):
            points.append({"id":i+1,"coords":(lcoords[i%4][0],lcoords[i%4][1],lcoords[i%4][2]-depth), "lc":None})


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

        # make line loops
        line_loops = []

        line_loops.append({"id":1,"lines":[1,2,3,4]})
        line_loops.append({"id":2,"lines":[-5,-8,-7,-6]})
        line_loops.append({"id":3,"lines":[10,6,-9,-2]})
        line_loops.append({"id":4,"lines":[-12,-4,11,8]})
        line_loops.append({"id":5,"lines":[-3,9,7,-11]})
        line_loops.append({"id":6,"lines":[5,-10,-1,12]})

        self.faces_names = {"up":5,"down":6,"right":3,"left":4,"front":1,"back":2}
        self.faces = {i: {"bubbles_in":[]} for i in range(1,7)}
        self.faces[1]["hyperplane"] = ((0,0,1), (0,0,self.max_z))
        self.faces[2]["hyperplane"] = ((0,0,-1),(0,0,self.min_z))
        self.faces[3]["hyperplane"] = ((1,0,0),(self.max_x,0,0))
        self.faces[4]["hyperplane"] = ((-1,0,0),(self.min_x,0,0))
        self.faces[5]["hyperplane"] = ((0,1,0),(0,self.max_y,0))
        self.faces[6]["hyperplane"] = ((0,-1,0),(0,self.min_y,0))

        
        self.geometry = {"points":{p["id"]: {"coords":p["coords"], "lc":None} for p in points}}
        self.geometry["lines"] = {l["id"]:(l["start"], l["end"]) for l in lines}
        self.geometry["lineLoops"] = {ll["id"]:ll["lines"] for ll in line_loops}



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
                data = bubble.discretize(accuracy, self.freePointId, self.freeLineId, self.freeLineLoopId)
                
                self.geometry["points"] = {**self.geometry["points"],**data["points"]}
                self.geometry["lines"] = {**self.geometry["lines"],**data["lines"]}
                lineLoops = []
                for triangle in data["lineLoops"]:

                    lineLoops.append({"id":triangle, "lines":data["lineLoops"][triangle]})

                self.geometry["lineLoops"] = {**self.geometry["lineLoops"], **data["lineLoops"]}
                self.interior_lineLoops += lineLoops

                self.freePointId += len(data["points"]) # TODO: correct this
                self.freeLineId += len(data["lines"])
                self.freeLineLoopId += len(lineLoops)

            elif loc in ["up","down","right","left","front","back"]:
                clippingPlane = self.faces[self.faces_names[loc]]["hyperplane"]
                cut_height = self.height_up
                bubble_geo = bubble.discretize(accuracy, self.freePointId, self.freeLineId, self.freeLineLoopId)


                clip_bubble, intersectionLoop = clip(bubble_geo, clippingPlane, EuclideanProjection, EasyClip)
                
                
                self.geometry["points"] = {**self.geometry["points"], **clip_bubble["points"]}
                self.geometry["lines"] = {**self.geometry["lines"], **clip_bubble["lines"]}
                self.geometry["lineLoops"] = {**self.geometry["lineLoops"],**clip_bubble["lineLoops"],**{intersectionLoop[0]:intersectionLoop[1]}}

                self.faces[self.faces_names[loc]]["bubbles_in"].append({"clip_geo":clip_bubble,"intersectionLoop":intersectionLoop})
                
                self.freePointId = (max([p for p in clip_bubble["points"]])+1)
                self.freeLineId = (max([l for l in clip_bubble["lines"]])+1)
                self.freeLineLoopId  = (max([ll for ll in clip_bubble["lineLoops"]])+2)
            else:
                raise NotImplementedError(loc)
        else:
            print("Bubble could not be added", loc, valid)
            
    def _finish_geometry(self):
        self.geometry["surfaces"] = dict()
        # Write exterior surfaces
        for face in self.faces:
            lintersectionLoops = []
            
            for bubble in self.faces[face]["bubbles_in"]:
                lintersectionLoops.append(bubble["intersectionLoop"][0])

            main_surface = {"id":face,"lineLoops":[face] + lintersectionLoops}
            
            
            self.exterior_surfaces.append(main_surface)
            
            self.geometry["surfaces"][face] = [face]+lintersectionLoops


            for bubble in self.faces[face]["bubbles_in"]:
                llineLoops = bubble["clip_geo"]["lineLoops"]
                
                for lineLoop in llineLoops:
                    Id = lineLoop
                    surface = {"id":Id,"lineLoops":[Id]}

                    self.geometry["surfaces"][lineLoop] = [lineLoop]
                    self.exterior_surfaces.append(surface)
            
        
        # Write interior surfaces
        for lineLoop in self.interior_lineLoops:
            self.geometry["surfaces"][lineLoop["id"]] = [lineLoop["id"]]
        # Write exterior surface loop
        self.geometry["surfaceLoops"] = {1:[surface["id"] for surface in self.exterior_surfaces]}
        
        # Write interior surface loop
        self.geometry["surfaceLoops"][2] = [lineLoop["id"] for lineLoop in self.interior_lineLoops]
        # Write volume
        self.geometry["volumes"] = {1:[1,2]}

    def to_geo(self, file_name):
        self._finish_geometry()
        write_geo(self.geometry, file_name)