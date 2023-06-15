import json
import random

from bubbles.gmsh_api import write_geo
from bubbles.three_d.bubble import Sphere, StarBubble, clip
from bubbles.three_d.clip import EasyClip, EuclideanProjection
from bubbles.utils import geometry_to_gmsh_entities


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
        self.geometry = {
            "points": dict(),
            "lines": dict(),
            "lineLoops": dict(),
            "surfaces": dict(),
            "surfaceLoops": {1: [], 2: []},  # 1 is outside, 0 is inside
            "volumes": {1: [1, 2]},
        }
        self.merge = {
            # here we assume that glue will be merged only on upper and lower side
            5: {
                "mainLoopLayer": {1: [], 2: [], 3: [], 4: []},
                "subLoops": [],
                "mainLoopGlue": {1: [], 2: [], 3: [], 4: []},
                "loopId": None,
                "surfaceId": None,
            },
            6: {
                "mainLoopLayer": {1: [], 2: [], 3: [], 4: []},
                "subLoops": [],
                "mainLoopGlue": {1: [], 2: [], 3: [], 4: []},
                "loopId": None,
                "surfaceId": None,
            },
        }

        self.max_x = lcoords[1][0]
        self.max_y = lcoords[2][1]
        self.max_z = 0
        self.min_x = lcoords[0][0]
        self.min_y = lcoords[0][1]
        self.min_z = -depth

        self.freePointId = 9
        self.freeLineId = 13
        self.freeLineLoopId = 1
        self.freeSurfaceId = 1

        # make points
        points = [{"id": i + 1, "coords": lcoords[i], "lc": None} for i in range(4)]

        for i in range(4, 8):
            points.append(
                {
                    "id": i + 1,
                    "coords": (
                        lcoords[i % 4][0],
                        lcoords[i % 4][1],
                        lcoords[i % 4][2] - depth,
                    ),
                    "lc": None,
                }
            )

        # make lines
        lines = []
        for i in range(4):
            lines.append({"id": i + 1, "start": i + 1, "end": (i + 1) % 4 + 1})

        for i in range(4, 8):
            lines.append({"id": i + 1, "start": i + 1, "end": (i + 1) % 4 + 4 + 1})

        lines.append({"id": 9, "start": 3, "end": 7})
        lines.append({"id": 10, "start": 2, "end": 6})

        lines.append({"id": 11, "start": 4, "end": 8})
        lines.append({"id": 12, "start": 1, "end": 5})

        self.faces = {i: {"innerLoops": []} for i in range(1, 7)}

        self.faces[1]["outerLoop"] = {
            1: {"adjacentFace": 6, "lines": [1], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 3, "lines": [2], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 5, "lines": [3], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 4, "lines": [4], "lineLoops": [], "connections": []},
        }
        self.faces[2]["outerLoop"] = {
            1: {"adjacentFace": 6, "lines": [-5], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 4, "lines": [-8], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 5, "lines": [-7], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 3, "lines": [-6], "lineLoops": [], "connections": []},
        }
        self.faces[3]["outerLoop"] = {
            1: {"adjacentFace": 6, "lines": [10], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 2, "lines": [6], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 5, "lines": [-9], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 1, "lines": [-2], "lineLoops": [], "connections": []},
        }
        self.faces[4]["outerLoop"] = {
            1: {"adjacentFace": 6, "lines": [-12], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 1, "lines": [-4], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 5, "lines": [11], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 2, "lines": [8], "lineLoops": [], "connections": []},
        }
        self.faces[5]["outerLoop"] = {
            1: {"adjacentFace": 1, "lines": [-3], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 3, "lines": [9], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 2, "lines": [7], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 4, "lines": [-11], "lineLoops": [], "connections": []},
        }
        self.faces[6]["outerLoop"] = {
            1: {"adjacentFace": 2, "lines": [5], "lineLoops": [], "connections": []},
            2: {"adjacentFace": 3, "lines": [-10], "lineLoops": [], "connections": []},
            3: {"adjacentFace": 1, "lines": [-1], "lineLoops": [], "connections": []},
            4: {"adjacentFace": 4, "lines": [12], "lineLoops": [], "connections": []},
        }
        for i in range(1, 7):
            for j in range(1, 5):
                self.faces[i]["outerLoop"][j]["edgeLines"] = {}

        self.faces[1]["hyperplane"] = ((0, 0, 1), (0, 0, self.max_z))
        self.faces[2]["hyperplane"] = ((0, 0, -1), (0, 0, self.min_z))
        self.faces[3]["hyperplane"] = ((1, 0, 0), (self.max_x, 0, 0))
        self.faces[4]["hyperplane"] = ((-1, 0, 0), (self.min_x, 0, 0))
        self.faces[5]["hyperplane"] = ((0, 1, 0), (0, self.max_y, 0))
        self.faces[6]["hyperplane"] = ((0, -1, 0), (0, self.min_y, 0))

        self.geometry["points"] = {
            p["id"]: {"coords": p["coords"], "lc": None} for p in points
        }
        self.geometry["lines"] = {l["id"]: (l["start"], l["end"]) for l in lines}
        # self.geometry["lineLoops"] = {ll:self.faces[ll]["outerLoop"] for ll in self.faces}
        # print(self.geometry["lineLoops"])
        self.ready = False  # Ready to make GeoFile

        self.epsilon = (
            min(
                self.max_x - self.min_x,
                self.max_y - self.min_y,
                self.max_z - self.min_z,
            )
            / 15
        )

    def sampleBubble(self):
        radius = random.uniform(
            0,
            min(
                self.max_x - self.min_x,
                self.max_y - self.min_y,
                self.max_z - self.min_z,
            )
            / 2,
        )
        mx = random.uniform(self.min_x + self.epsilon, self.max_x - self.epsilon)
        my = random.uniform(self.min_y + self.epsilon, self.max_y - self.epsilon)
        mz = random.uniform(self.min_z + self.epsilon, self.max_z - self.epsilon)
        s = Sphere(radius, (mx, my, mz))
        self.add_bubble(s)

    def add_bubble(self, bubble, accuracy=0.5):
        loc, valid = self._localize_bubble(bubble)
        if valid:
            if len(loc) == 0:
                print("Bubble in the interior of the glue")
                # Bubble inside
                data = bubble.discretize(
                    accuracy, self.freePointId, self.freeLineId, self.freeLineLoopId
                )

                self.geometry["points"] = {**self.geometry["points"], **data["points"]}
                self.geometry["lines"] = {**self.geometry["lines"], **data["lines"]}
                lineLoops = []
                for triangle in data["lineLoops"]:
                    lineLoops.append(
                        {"id": triangle, "lines": data["lineLoops"][triangle]}
                    )
                    self.geometry["lineLoops"][triangle] = data["lineLoops"][triangle]
                    self.geometry["surfaces"][triangle] = [triangle]
                    self.geometry["surfaceLoops"][2].append(triangle)

                self.geometry["lineLoops"] = {
                    **self.geometry["lineLoops"],
                    **data["lineLoops"],
                }

                self.freePointId += len(data["points"])  # TODO: correct this
                self.freeLineId += len(data["lines"])
                self.freeLineLoopId += len(lineLoops)

            elif len(loc) == 1:
                print("Bubble intersects with a facet!")
                # Intersection with face
                loc = loc[0]
                clippingPlane = self.faces[loc]["hyperplane"]
                bubble_geo = bubble.discretize(
                    accuracy, self.freePointId, self.freeLineId, self.freeLineLoopId
                )

                clip_bubble, intersectionLoop = clip(
                    bubble_geo, clippingPlane, EuclideanProjection, EasyClip
                )

                self.geometry["points"] = {
                    **self.geometry["points"],
                    **clip_bubble["points"],
                }
                self.geometry["lines"] = {
                    **self.geometry["lines"],
                    **clip_bubble["lines"],
                }
                self.geometry["lineLoops"] = {
                    **self.geometry["lineLoops"],
                    **clip_bubble["lineLoops"],
                    **{intersectionLoop[0]: intersectionLoop[1]},
                }

                for triangle in clip_bubble["lineLoops"]:
                    self.geometry["surfaces"][triangle] = [triangle]
                    self.geometry["surfaceLoops"][1].append(triangle)

                self.faces[loc]["innerLoops"].append(intersectionLoop[0])
                if loc in self.merge:
                    self.merge[loc]["subLoops"].append(intersectionLoop[0])
                self.freePointId = max([p for p in clip_bubble["points"]]) + 1
                self.freeLineId = max([l for l in clip_bubble["lines"]]) + 1
                self.freeLineLoopId = (
                    intersectionLoop[0] + 1
                )  # (max([ll for ll in clip_bubble["lineLoops"]])+2)
                self.freeSurfaceId = self.freeLineLoopId

            elif len(loc) == 2:
                print("Bubble intersects with an edge!")
                # Intersection with edge
                clippingPlane0, clippingPlane1 = (
                    self.faces[loc[0]]["hyperplane"],
                    self.faces[loc[1]]["hyperplane"],
                )
                bubble_geo = bubble.discretize(
                    accuracy, self.freePointId, self.freeLineId, self.freeLineLoopId
                )
                clip_bubble0, intersectionLoop0 = clip(
                    bubble_geo, clippingPlane0, EuclideanProjection, EasyClip
                )
                clip_bubble1, intersectionLoop1 = clip(
                    clip_bubble0, clippingPlane1, EuclideanProjection, EasyClip
                )

                loop0_clip0 = []
                loop0_clip1 = []
                out = False
                for i, l in enumerate(intersectionLoop0[1]):
                    if abs(l) in clip_bubble1["lines"] and not out:
                        loop0_clip0.append(l)
                    elif abs(l) in clip_bubble1["lines"] and out:
                        loop0_clip1.append(l)
                    else:
                        out = True
                loop0 = loop0_clip1 + loop0_clip0
                intersectionLoop0 = intersectionLoop0[0], loop0
                intersectionLoops = (intersectionLoop0, intersectionLoop1)

                # add loop to face, find correct side
                for i in range(2):
                    localSide = None
                    for side, val in self.faces[loc[i]]["outerLoop"].items():
                        if val["adjacentFace"] == loc[1 - i]:
                            localSide = side
                    # Check if the orientation of intersectionLoop is correct
                    line_id = abs(
                        self.faces[loc[i]]["outerLoop"][localSide]["lines"][0]
                    )
                    lineLoop = intersectionLoops[i][1]
                    startPoint_Loop = (
                        clip_bubble1["lines"][lineLoop[0]][0]
                        if lineLoop[0] > 0
                        else clip_bubble1["lines"][-lineLoop[0]][1]
                    )

                    endPoint_Loop = (
                        clip_bubble1["lines"][lineLoop[-1]][1]
                        if lineLoop[-1] > 0
                        else clip_bubble1["lines"][-lineLoop[-1]][0]
                    )
                    point0 = clip_bubble1["points"][startPoint_Loop]["coords"]
                    point1 = clip_bubble1["points"][endPoint_Loop]["coords"]
                    l0, l1 = self.geometry["lines"][line_id]
                    l0, l1 = (
                        self.geometry["points"][l0]["coords"],
                        self.geometry["points"][l1]["coords"],
                    )
                    if self.faces[loc[i]]["outerLoop"][localSide]["lines"][0] < 0:
                        temp = l1
                        l1 = l0
                        l0 = temp
                    if not compare(point0, point1, (l0, l1)):
                        self.faces[loc[i]]["outerLoop"][localSide]["lineLoops"].append(
                            [-l for l in reversed(intersectionLoops[i][1])]
                        )
                    else:
                        self.faces[loc[i]]["outerLoop"][localSide]["lineLoops"].append(
                            intersectionLoops[i][1]
                        )

                self.geometry["points"] = {
                    **self.geometry["points"],
                    **clip_bubble1["points"],
                }
                self.geometry["lines"] = {
                    **self.geometry["lines"],
                    **clip_bubble1["lines"],
                }
                self.geometry["lineLoops"] = {
                    **self.geometry["lineLoops"],
                    **clip_bubble1["lineLoops"],
                }

                for triangle in clip_bubble1["lineLoops"]:
                    self.geometry["surfaces"][triangle] = [triangle]
                    self.geometry["surfaceLoops"][1].append(triangle)

                self.freePointId = max([p for p in clip_bubble1["points"]]) + 1
                self.freeLineId = max([l for l in clip_bubble1["lines"]]) + 1
                self.freeLineLoopId = intersectionLoop1[0] + 1
                self.freeSurfaceId = self.freeLineLoopId
            else:
                raise NotImplementedError(loc)
        else:
            print("Bubble could not be added", loc, valid)

    def _finish_geometry(self):
        """
        After Bubbles have been added to the Box write LineLoops defining the 6 faces.
        """
        # Write exterior surface loop

        JunkLines = set()  # lines to be removed
        # Iterate over faces
        for i in range(1, 7):
            faceLineLoop = []
            # Iterate over Edges
            for j in range(1, 5):
                edgeLines = []
                edgeLoop = []
                if (
                    len(self.faces[i]["outerLoop"][j]["lineLoops"]) == 0
                ):  # No intersection with Edges
                    lineID = self.faces[i]["outerLoop"][j]["lines"][0]
                    edgeLoop.append(lineID)
                    edgeLines.append((lineID, self.geometry["lines"][abs(lineID)]))
                else:
                    line_id = abs(self.faces[i]["outerLoop"][j]["lines"][0])
                    JunkLines.add(line_id)
                    # sort list of LineLoops for correct orientation, anticlockwise
                    self.faces[i]["outerLoop"][j]["lineLoops"].sort(
                        key=sortLineLoops(
                            self.faces[i]["outerLoop"][j]["lines"][0], self.geometry
                        )
                    )

                    if (
                        len(self.faces[i]["outerLoop"][j]["connections"]) == 0
                    ):  # len(self.faces[i]["outerLoop"][j]["connections"]) == 0:
                        # write connection Lines, e.g. lines which connect the holes
                        connections = []

                        startEndPoints = []
                        firstLine = self.faces[i]["outerLoop"][j]["lines"][0]
                        firstPoint = (
                            self.geometry["lines"][firstLine][0]
                            if firstLine > 0
                            else self.geometry["lines"][-firstLine][1]
                        )
                        for lineLoop in self.faces[i]["outerLoop"][j]["lineLoops"]:
                            firstLine, lastLine = lineLoop[0], lineLoop[-1]
                            lastPoint = (
                                self.geometry["lines"][firstLine][0]
                                if firstLine > 0
                                else self.geometry["lines"][-firstLine][1]
                            )
                            startEndPoints.append((firstPoint, lastPoint))
                            firstPoint = (
                                self.geometry["lines"][lastLine][1]
                                if lastLine > 0
                                else self.geometry["lines"][-lastLine][0]
                            )
                        lastLine = self.faces[i]["outerLoop"][j]["lines"][0]
                        lastPoint = (
                            self.geometry["lines"][lastLine][1]
                            if lastLine > 0
                            else self.geometry["lines"][-lastLine][0]
                        )
                        startEndPoints.append((firstPoint, lastPoint))
                        # Write new lines and add all lines to faceLineLoop
                        assert (
                            len(startEndPoints)
                            == len(self.faces[i]["outerLoop"][j]["lineLoops"]) + 1
                        )
                        for m in range(len(startEndPoints) - 1):
                            newLine = startEndPoints[m]
                            newLineID = self.freeLineId
                            self.geometry["lines"][newLineID] = newLine
                            connections.append(newLineID)
                            self.freeLineId += 1
                            #
                            edgeLoop.append(newLineID)

                            edgeLoop = (
                                edgeLoop + self.faces[i]["outerLoop"][j]["lineLoops"][m]
                            )
                            edgeLines.append(
                                (newLineID, self.geometry["lines"][abs(newLineID)])
                            )
                            a = self.faces[i]["outerLoop"][j]["lineLoops"][m][0]
                            a = (
                                self.geometry["lines"][a][0]
                                if a > 0
                                else self.geometry["lines"][-a][1]
                            )
                            b = self.faces[i]["outerLoop"][j]["lineLoops"][m][-1]
                            b = (
                                self.geometry["lines"][b][1]
                                if b > 0
                                else self.geometry["lines"][-b][0]
                            )
                            newLineMerge = (a, b)
                            edgeLines.append((None, (a, b)))
                        newLine = startEndPoints[-1]
                        newLineID = self.freeLineId
                        self.geometry["lines"][newLineID] = newLine
                        connections.append(newLineID)
                        self.freeLineId += 1

                        edgeLoop.append(newLineID)
                        edgeLines.append(
                            (newLineID, self.geometry["lines"][abs(newLineID)])
                        )

                        adjacentFaceId = self.faces[i]["outerLoop"][j]["adjacentFace"]
                        t = None
                        for key, val in self.faces[adjacentFaceId]["outerLoop"].items():
                            if val["adjacentFace"] == i:
                                t = key
                        self.faces[adjacentFaceId]["outerLoop"][t]["connections"] = [
                            -l for l in reversed(connections)
                        ]

                    else:
                        adjacentFaceId = self.faces[i]["outerLoop"][j]["adjacentFace"]
                        t = None
                        for key, val in self.faces[adjacentFaceId]["outerLoop"].items():
                            if val["adjacentFace"] == i:
                                t = key
                        self.faces[adjacentFaceId]["outerLoop"][t]["connections"] = [
                            -l for l in reversed(connections)
                        ]
                        connections = self.faces[adjacentFaceId]["outerLoop"][t][
                            "connections"
                        ]
                        for m in range(len(connections) - 1):
                            edgeLoop.append(connections[m])

                            edgeLoop = (
                                edgeLoop + self.faces[i]["outerLoop"][j]["lineLoops"][m]
                            )
                            edgeLines.append(
                                (
                                    connections[m],
                                    self.geometry["lines"][abs(connections[m])],
                                )
                            )

                            a = self.faces[i]["outerLoop"][j]["lineLoops"][m][0]
                            a = (
                                self.geometry["lines"][a][0]
                                if a > 0
                                else self.geometry["lines"][-a][1]
                            )
                            b = self.faces[i]["outerLoop"][j]["lineLoops"][m][-1]
                            b = (
                                self.geometry["lines"][b][1]
                                if b > 0
                                else self.geometry["lines"][-b][0]
                            )
                            newLineMerge = (a, b)
                            edgeLines.append((None, (a, b)))

                        edgeLoop.append(connections[-1])
                        edgeLines.append(
                            (
                                connections[-1],
                                self.geometry["lines"][abs(connections[-1])],
                            )
                        )

                self.faces[i]["outerLoop"][j]["edgeLines"] = edgeLines
                faceLineLoop = faceLineLoop + edgeLoop
                if i in {5, 6}:
                    adjacentFaceId = self.faces[i]["outerLoop"][j]["adjacentFace"]
                    self.merge[i]["mainLoopLayer"][adjacentFaceId] = edgeLines
                    self.merge[i]["mainLoopGlue"][adjacentFaceId] = edgeLoop

            self.geometry["lineLoops"][self.freeLineLoopId] = faceLineLoop
            self.geometry["surfaces"][self.freeSurfaceId] = [
                self.freeLineLoopId
            ] + self.faces[i]["innerLoops"]
            if i in {5, 6}:
                self.merge[i]["loopId"] = self.freeLineLoopId
                self.merge[i]["surfaceId"] = self.freeSurfaceId
            self.freeLineLoopId += 1
            self.geometry["surfaceLoops"][1].append(self.freeSurfaceId)
            self.freeSurfaceId += 1

        for lineId in JunkLines:
            self.geometry["lines"].pop(lineId)
        self.ready = True

    def _localize_bubble(self, bubble):
        """
        front: 1
        back: 2
        right:3
        left:4
        up:5
        down:6
        """

        # Assume midpoint is inside of the box and the bubble is small enough not to intersect with opposite facets

        max_x, max_y, max_z, min_x, min_y, min_z = bubble.getBoundingBox()
        valid = True
        loc = []
        if max_x >= self.max_x:
            # intersection with right
            loc.append(3)
            if max_y >= self.max_y:
                # intersection with up
                loc.append(5)
                if max_z >= self.max_z:
                    # intersection with front
                    # cornercase
                    loc.append(1)
                    valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    # cornercase
                    loc.append(2)
                    valid = False
                else:
                    pass  # valid = False
                # no intersection with right/left
            elif min_y <= self.min_y:
                # intersection with down
                loc.append(6)
                if max_z >= self.max_z:
                    # intersection with front
                    # cornercase
                    loc.append(1)
                    valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    # cornercase
                    loc.append(2)
                    valid = False
                else:
                    # no intersection with front/back
                    pass  # valid = False
            else:
                # no intersection with up/down
                if max_z >= self.max_z:
                    # intersection with front
                    loc.append(1)
                    pass  # valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    loc.append(2)
                    pass  # valid = False
                else:
                    # no intersection with front/back
                    valid = True

        elif min_x <= self.min_x:
            # intersection with left
            loc.append(4)
            if max_y >= self.max_y:
                # intersection with up
                loc.append(5)
                if max_z >= self.max_z:
                    # intersection with front
                    # cornercase
                    loc.append(1)
                    valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    # cornercase
                    loc.append(2)
                    valid = False
                else:
                    # no intersection with front/back
                    # valid = False
                    pass
            elif min_y <= self.min_y:
                # intersection with down
                loc.append(6)
                if max_z >= self.max_z:
                    # intersection with front
                    # cornercase
                    loc.append(1)
                    valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    # cornercase
                    loc.append(2)
                    valid = False
                else:
                    # no intersection with front/back
                    # valid = False
                    pass

            else:
                # no intersection with up/down
                if max_z >= self.max_z:
                    # intersection with front
                    loc.append(1)
                    pass  # valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    loc.append(2)
                    pass  # valid = False
                else:
                    # no intersection with front/back
                    pass

        else:
            # no intersection with left/right
            if max_y >= self.max_y:
                # intersection with up
                loc.append(5)
                if max_z >= self.max_z:
                    # intersection with front
                    loc.append(1)
                    pass  # valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    loc.append(2)
                    pass  # valid = False
                else:
                    # no intersection with front/back
                    pass
            elif min_y <= self.min_y:
                # intersection with down
                loc.append(6)
                if max_z >= self.max_z:
                    # intersection with front
                    loc.append(1)
                    pass  # valid = False
                elif min_z <= self.min_z:
                    # intersection with back
                    loc.append(2)
                    pass  # valid = False
                else:
                    # no intersection with front/back
                    pass
            else:
                # no intersection with up/down
                if max_z >= self.max_z:
                    # intersection with front
                    loc.append(1)
                elif min_z <= self.min_z:
                    # intersection with back
                    loc.append(2)
                else:
                    # no intersection with anyside
                    pass

        return loc, valid

    def to_geo(self, file_name):
        self._finish_geometry()
        gmsh_entity_dict = geometry_to_gmsh_entities(self.geometry)
        write_geo(file_name, **gmsh_entity_dict)

    def showFaces(self):
        print(json.dumps(self.faces, indent=4))


def compare(point0, point1, line):
    start, end = line
    v = [end[i] - start[i] for i in range(3)]
    p0 = [point0[i] - start[i] for i in range(3)]
    p1 = [point1[i] - start[i] for i in range(3)]

    val0 = sum([p0[i] * v[i] for i in range(3)])
    val1 = sum([p1[i] * v[i] for i in range(3)])
    return val0 <= val1


def sortLineLoops(line, geometry):
    start, end = (
        geometry["lines"][line]
        if line > 0
        else (geometry["lines"][-line][1], geometry["lines"][-line][0])
    )
    startPoint = geometry["points"][start]["coords"]
    endPoint = geometry["points"][end]["coords"]

    def f(ll):
        firstLine = ll[0]
        firstPoint = (
            geometry["lines"][firstLine][0]
            if firstLine > 0
            else geometry["lines"][-firstLine][1]
        )
        firstPoint = geometry["points"][firstPoint]["coords"]
        v = [endPoint[i] - startPoint[i] for i in range(3)]
        p0 = [firstPoint[i] - startPoint[i] for i in range(3)]
        val0 = sum([p0[i] * v[i] for i in range(3)])

        return val0

    return f
