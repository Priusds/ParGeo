from bubbles.geo_utils.utils import write_geo
from bubbles.three_d.glue import *
from bubbles.three_d.layers import *


def get_key(val, my_dict):
    for key, value in my_dict.items():
        if val == value:
            return key

    raise ValueError("No key in my_dict with value: ", val)


def replacePoint(Id, newPoint, geometry):
    newID = newPoint[0]
    pointDict = newPoint[1]

    geometry["points"].pop(Id, None)
    geometry["points"][newID] = pointDict
    for lineID, points in geometry["lines"].items():
        if Id == points[0]:
            geometry["lines"][lineID] = (newID, points[1])
        elif Id == points[1]:
            geometry["lines"][lineID] = (points[0], newID)
        else:
            pass
    return geometry


def replaceLine(IdOld, newLines, geometry):
    # print(geometry)
    assert IdOld > 0
    oldLine = geometry["lines"][IdOld]

    if len(newLines) > 2:
        lnewID = newLines

        lnewlines = [
            geometry["lines"][Id] if Id > 0 else geometry["lines"][-Id]
            for Id in newLines
        ]
    else:
        lnewID = list(newLines[0])

        lnewlines = list(newLines[1])

    # Check if new and old line have same orientation
    startOld = oldLine[0]
    startNew = lnewlines[0][0] if lnewID[0] > 0 else lnewlines[0][1]
    sameOrientation = startOld == startNew

    if not sameOrientation:
        lnewID = [-l if l is not None else None for l in reversed(lnewID)]
        lnewlines.reverse()

    geometry["lines"].pop(IdOld, None)

    for i, newID in enumerate(lnewID):  # add new lines
        if newID is None:
            freeLineID = max(geometry["lines"].keys()) + 1

            lnewID[i] = freeLineID

            endPoint = (
                geometry["lines"][lnewID[i - 1]][1]
                if lnewID[i - 1] > 0
                else geometry["lines"][-lnewID[i - 1]][0]
            )
            if endPoint != lnewlines[i][0]:
                lnewlines[i] = (endPoint, lnewlines[i][0])
            geometry["lines"][freeLineID] = lnewlines[i]
        else:
            geometry["lines"][abs(newID)] = lnewlines[i]

    for ll, llines in geometry["lineLoops"].items():
        if IdOld in llines:
            idx = llines.index(IdOld)
            geometry["lineLoops"][ll] = llines[0:idx] + lnewID + llines[idx + 1 :]

        elif -IdOld in llines:
            idx = llines.index(-IdOld)
            geometry["lineLoops"][ll] = (
                llines[0:idx] + [-l for l in reversed(lnewID)] + llines[idx + 1 :]
            )

        else:
            pass

    return lnewID if sameOrientation else [-l for l in reversed(lnewID)]


def replaceLineLoop(oldID, newLineLoop, geometry):
    # TODO: check if same orientation, if not add -* to surfaces.
    if type(newLineLoop) == int:
        newID = newLineLoop
        newLines = geometry["lineLoops"][newID]
    elif len(newLineLoop) == 2:
        newID = newLineLoop[0]
        newLines = newLineLoop[1]

    assert newID > 0

    geometry["lineLoops"].pop(oldID)
    geometry["lineLoops"][newID] = newLines
    for surfaceID, lineLoops in geometry["surfaces"].items():
        if oldID in lineLoops:
            idx = lineLoops.index(oldID)
            geometry["surfaces"][surfaceID] = (
                lineLoops[0:idx] + [newID] + lineLoops[idx + 1 :]
            )
        elif -oldID in lineLoops:
            idx = lineLoops.index(-oldID)
            geometry["surfaces"][surfaceID] = (
                lineLoops[0:idx] + [-newID] + lineLoops[idx + 1 :]
            )
        else:
            pass


def replaceSurface(oldID, newSurface, geometry):
    if type(newSurface) == int:
        newID = newSurface
        newLineLoops = geometry["surfaces"][newID]
    elif len(newLineLoop) == 2:
        newID = newSurface[0]
        newLineLoops = newSurface[1]

    assert newID > 0
    # assert oldID in geometry["surfaces"]

    geometry["surfaces"].pop(oldID)
    geometry["surfaces"][newID] = newLineLoops
    for sLoopID, surfaces in geometry["surfaceLoops"].items():
        if oldID in surfaces:
            idx = surfaces.index(oldID)
            geometry["surfaceLoops"][sLoopID] = (
                surfaces[0:idx] + [newID] + surfaces[idx + 1 :]
            )
        elif -oldID in surfaces:
            idx = surfaces.index(-oldID)
            geometry["surfaceLoops"][sLoopID] = (
                surfaces[0:idx] + [-newID] + surfaces[idx + 1 :]
            )
        else:
            pass


class Sandwich:
    def __init__(self, layers1, glue, layers2):
        """
        444 LAYER2 333

         GLUE

        44 LAYER1 33

        Assume we see the front side of the sandwich
        """
        assert glue.ready == True
        self.nLayers = (
            layers1.nLayers
        )  # Assume both layers are made of the same number of sublayers
        glueGeometry = glue.geometry
        layers1Geometry = layers1.geometry
        layers2Geometry = layers2.geometry

        self.merge = {
            33: layers1.merge[3],
            44: layers1.merge[4],
            333: layers2.merge[3],
            444: layers2.merge[4],
        }

        sandwich = {
            "points": glueGeometry["points"],
            "lines": glueGeometry["lines"],
            "lineLoops": glueGeometry["lineLoops"],
            "surfaces": glueGeometry["surfaces"],
            "surfaceLoops": glueGeometry["surfaceLoops"],
            "volumes": glueGeometry["volumes"],
        }

        maxGlueId = max(
            max(glueGeometry["points"].keys()),
            max(glueGeometry["lines"].keys()),
            max(glueGeometry["lineLoops"].keys()),
            max(glueGeometry["surfaces"].keys()),
        )

        maxLayer2Id = max(
            max(layers2Geometry["points"].keys()),
            max(layers2Geometry["lines"].keys()),
            max(layers2Geometry["lineLoops"].keys()),
            max(layers2Geometry["surfaces"].keys()),
        )
        # ============================================
        # NOW SHIFT LAYERS IDS for the merge
        # ============================================

        # Layer2 gets shifted by maxGlueId
        # Layer1 gets shifted by maxGlueId + maxLayer2Id

        for i in range(2):
            shiftFactor = [maxGlueId, maxGlueId + maxLayer2Id][i]
            layer = [layers2Geometry, layers1Geometry][i]

            layer["points"] = {
                shiftFactor + Id: layer["points"][Id] for Id in layer["points"].keys()
            }
            layer["lines"] = {
                shiftFactor
                + Id: tuple(pointId + shiftFactor for pointId in layer["lines"][Id])
                for Id in layer["lines"].keys()
            }
            layer["lineLoops"] = {
                shiftFactor
                + Id: [
                    lineId + shiftFactor if lineId > 0 else lineId - shiftFactor
                    for lineId in layer["lineLoops"][Id]
                ]
                for Id in layer["lineLoops"].keys()
            }
            layer["surfaces"] = {
                shiftFactor
                + Id: [lineLoopId + shiftFactor for lineLoopId in layer["surfaces"][Id]]
                for Id in layer["surfaces"].keys()
            }
            layer["surfaceLoops"] = {
                shiftFactor
                + Id: [
                    surfaceId + shiftFactor for surfaceId in layer["surfaceLoops"][Id]
                ]
                for Id in layer["surfaceLoops"].keys()
            }
            layer["volumes"] = {
                shiftFactor
                + Id: [
                    surfaceLoopId + shiftFactor
                    for surfaceLoopId in layer["volumes"][Id]
                ]
                for Id in layer["volumes"].keys()
            }
            sandwich = {
                "points": {**sandwich["points"], **layer["points"]},
                "lines": {**sandwich["lines"], **layer["lines"]},
                "lineLoops": {**sandwich["lineLoops"], **layer["lineLoops"]},
                "surfaces": {**sandwich["surfaces"], **layer["surfaces"]},
                "surfaceLoops": {**sandwich["surfaceLoops"], **layer["surfaceLoops"]},
                "volumes": {**sandwich["volumes"], **layer["volumes"]},
            }

            self.merge[[444, 44][i]]["points"]["left"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["points"]["left"]
            ]
            self.merge[[333, 33][i]]["points"]["left"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["points"]["left"]
            ]
            self.merge[[444, 44][i]]["points"]["right"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["points"]["right"]
            ]
            self.merge[[333, 33][i]]["points"]["right"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["points"]["right"]
            ]

            self.merge[[444, 44][i]]["lines"]["left"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["lines"]["left"]
            ]
            self.merge[[333, 33][i]]["lines"]["left"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["lines"]["left"]
            ]
            self.merge[[444, 44][i]]["lines"]["right"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["lines"]["right"]
            ]
            self.merge[[333, 33][i]]["lines"]["right"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["lines"]["right"]
            ]
            self.merge[[444, 44][i]]["lines"]["middle"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["lines"]["middle"]
            ]
            self.merge[[333, 33][i]]["lines"]["middle"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["lines"]["middle"]
            ]

            self.merge[[333, 33][i]]["lineLoops"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["lineLoops"]
            ]
            self.merge[[444, 44][i]]["lineLoops"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["lineLoops"]
            ]

            self.merge[[333, 33][i]]["surfaces"] = [
                Id + shiftFactor for Id in self.merge[[333, 33][i]]["surfaces"]
            ]
            self.merge[[444, 44][i]]["surfaces"] = [
                Id + shiftFactor for Id in self.merge[[444, 44][i]]["surfaces"]
            ]

        freeLineId = max(sandwich["lines"].keys()) + 1
        freeLineLoopId = max(sandwich["lineLoops"].keys()) + 1
        freeSurfaceId = max(sandwich["surfaces"].keys()) + 1
        freeSurfaceLoopId = max(sandwich["surfaceLoops"].keys()) + 1
        freeVolumeId = max(sandwich["volumes"].keys()) + 1

        # =============================================
        # NOW ELIMINATE DUPLICATES
        # =============================================
        case = 2
        if case == 1:  # Glue and layers are aligned
            # POINTS
            gluePointIdsUp = [4, 3, 7, 8]
            layers2PointsId = [Id + maxGlueId for Id in layers2.merge[6]["points"]]

            gluePointIdsDown = [1, 2, 6, 5]
            layers1PointsId = [
                Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["points"]
            ]

            for i in range(4):
                replacePoint(
                    layers2PointsId[i],
                    (gluePointIdsUp[i], sandwich["points"][gluePointIdsUp[i]]),
                    sandwich,
                )
                replacePoint(
                    layers1PointsId[i],
                    (gluePointIdsDown[i], sandwich["points"][gluePointIdsDown[i]]),
                    sandwich,
                )

            self.merge[444]["points"]["right"][0] = 4
            self.merge[333]["points"]["left"][0] = 3
            self.merge[333]["points"]["right"][0] = 7
            self.merge[444]["points"]["left"][0] = 8

            self.merge[44]["points"]["right"][-1] = 1
            self.merge[33]["points"]["left"][-1] = 2
            self.merge[33]["points"]["right"][-1] = 6
            self.merge[44]["points"]["left"][-1] = 5

            # LINES
            layers2LinesId = [Id + maxGlueId for Id in layers2.merge[6]["lines"]]
            layers1LinesId = [
                Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["lines"]
            ]

            replaceLine(
                layers2LinesId[0],
                list(zip(*glue.merge[5]["mainLoopLayer"][1])),
                sandwich,
            )
            replaceLine(
                layers2LinesId[1],
                list(zip(*glue.merge[5]["mainLoopLayer"][3])),
                sandwich,
            )
            replaceLine(
                layers2LinesId[2],
                list(zip(*glue.merge[5]["mainLoopLayer"][2])),
                sandwich,
            )
            replaceLine(
                layers2LinesId[3],
                list(zip(*glue.merge[5]["mainLoopLayer"][4])),
                sandwich,
            )

            self.merge[333]["lines"]["middle"][0] = abs(
                glue.merge[5]["mainLoopLayer"][3][0][0]
            )
            self.merge[444]["lines"]["middle"][0] = abs(
                glue.merge[5]["mainLoopLayer"][4][0][0]
            )

            replaceLine(
                layers1LinesId[0],
                list(zip(*glue.merge[6]["mainLoopLayer"][1])),
                sandwich,
            )
            replaceLine(
                layers1LinesId[1],
                list(zip(*glue.merge[6]["mainLoopLayer"][3])),
                sandwich,
            )
            replaceLine(
                layers1LinesId[2],
                list(zip(*glue.merge[6]["mainLoopLayer"][2])),
                sandwich,
            )
            replaceLine(
                layers1LinesId[3],
                list(zip(*glue.merge[6]["mainLoopLayer"][4])),
                sandwich,
            )

            self.merge[33]["lines"]["middle"][-1] = abs(
                glue.merge[6]["mainLoopLayer"][3][0][0]
            )
            self.merge[44]["lines"]["middle"][-1] = abs(
                glue.merge[6]["mainLoopLayer"][4][0][0]
            )

            # LINELOOPS

            replaceLineLoop(
                layers2.merge[6]["lineLoop"] + maxGlueId,
                glue.merge[5]["loopId"],
                sandwich,
            )

            replaceLineLoop(
                layers1.merge[5]["lineLoop"] + maxGlueId + maxLayer2Id,
                glue.merge[6]["loopId"],
                sandwich,
            )

            # SURFACES
            replaceSurface(
                layers2.merge[6]["surface"] + maxGlueId,
                glue.merge[5]["surfaceId"],
                sandwich,
            )
            replaceSurface(
                layers1.merge[5]["surface"] + maxGlueId + maxLayer2Id,
                glue.merge[6]["surfaceId"],
                sandwich,
            )

            # Now we need to define new surfaces for the intersection between holes and facet
            freeSurfaceId = max(sandwich["surfaces"].keys()) + 1
            subSurfacesUp = list(
                range(freeSurfaceId, freeSurfaceId + len(glue.merge[5]["subLoops"]))
            )
            for sl in glue.merge[5]["subLoops"]:
                sandwich["surfaces"][freeSurfaceId] = [sl]
                freeSurfaceId += 1

            subSurfacesDown = list(
                range(freeSurfaceId, freeSurfaceId + len(glue.merge[6]["subLoops"]))
            )
            for sl in glue.merge[6]["subLoops"]:
                sandwich["surfaces"][freeSurfaceId] = [sl]
                freeSurfaceId += 1
            sandwich["surfaceLoops"][layers2.merge[6]["surfaceLoop"] + maxGlueId] = (
                sandwich["surfaceLoops"][layers2.merge[6]["surfaceLoop"] + maxGlueId]
                + subSurfacesUp
            )
            sandwich["surfaceLoops"][
                layers1.merge[5]["surfaceLoop"] + maxGlueId + maxLayer2Id
            ] = (
                sandwich["surfaceLoops"][
                    layers1.merge[5]["surfaceLoop"] + maxGlueId + maxLayer2Id
                ]
                + subSurfacesDown
            )

        elif case == 2:  # case Glue inside sandwich
            # Remove layer Lines

            lineLeftFront = freeLineId
            sandwich["lines"][freeLineId] = (
                4,
                [Id + maxGlueId for Id in layers2.merge[6]["points"]][0],
            )
            freeLineId += 1

            lineRightFront = freeLineId
            sandwich["lines"][freeLineId] = (
                [Id + maxGlueId for Id in layers2.merge[6]["points"]][1],
                3,
            )
            freeLineId += 1

            lineRightBack = freeLineId

            sandwich["lines"][freeLineId] = (
                [Id + maxGlueId for Id in layers2.merge[6]["points"]][2],
                7,
            )
            freeLineId += 1

            lineLeftBack = freeLineId
            sandwich["lines"][freeLineId] = (
                8,
                [Id + maxGlueId for Id in layers2.merge[6]["points"]][3],
            )
            freeLineId += 1

            replaceLine(
                [Id + maxGlueId for Id in layers2.merge[6]["lines"]][0],
                [-lineLeftFront, -3, -lineRightFront],
                sandwich,
            )
            replaceLine(
                [Id + maxGlueId for Id in layers2.merge[6]["lines"]][2],
                [lineRightBack, 7, lineLeftBack],
                sandwich,
            )

            # Remove lineLoop
            sandwich["lineLoops"].pop(maxGlueId + layers2.merge[6]["lineLoop"])

            # Remove surface
            sandwich["surfaces"].pop(maxGlueId + layers2.merge[6]["surface"])
            sandwich["surfaceLoops"][
                layers2.merge[6]["surfaceLoop"] + maxGlueId
            ].remove(maxGlueId + layers2.merge[6]["surface"])

            # make lineLoops

            lineLoopLeft = glue.merge[5]["mainLoopGlue"][4] + [
                lineLeftFront,
                [Id + maxGlueId for Id in layers2.merge[6]["lines"]][-1],
                -lineLeftBack,
            ]
            # print(glue.merge[5]["mainLoopGlue"])
            # print(list(list(zip(*glue.merge[5]["mainLoopLayer"][3]))[0]))
            lineLoopRight = glue.merge[5]["mainLoopGlue"][3] + [
                -lineRightBack,
                -[Id + maxGlueId for Id in layers2.merge[6]["lines"]][1],
                lineRightFront,
            ]
            sandwich["lineLoops"][freeLineLoopId] = lineLoopLeft
            freeLineLoopId += 1

            sandwich["lineLoops"][freeLineLoopId] = lineLoopRight
            freeLineLoopId += 1

            # make surfaces
            surfaceLeft = freeSurfaceId
            sandwich["surfaces"][surfaceLeft] = [freeLineLoopId - 2]
            freeSurfaceId += 1

            surfaceRight = freeSurfaceId
            sandwich["surfaces"][surfaceRight] = [freeLineLoopId - 1]
            freeSurfaceId += 1

            # Now we need to define new surfaces for the intersection between holes and facet
            subSurfacesUp = list(
                range(freeSurfaceId, freeSurfaceId + len(glue.merge[5]["subLoops"]))
            )
            for sl in glue.merge[5]["subLoops"]:
                sandwich["surfaces"][freeSurfaceId] = [sl]
                freeSurfaceId += 1

            # make SurfaceLoop
            sandwich["surfaceLoops"][layers2.merge[6]["surfaceLoop"] + maxGlueId] = (
                sandwich["surfaceLoops"][layers2.merge[6]["surfaceLoop"] + maxGlueId]
                + [surfaceLeft, surfaceRight]
                + subSurfacesUp
                + [glue.merge[5]["surfaceId"]]
            )

            # ========================================
            # LOWER PART
            # ========================================

            # Remove layer Lines

            lineLeftFront = freeLineId
            sandwich["lines"][freeLineId] = (
                1,
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["points"]][0],
            )
            freeLineId += 1

            lineRightFront = freeLineId
            sandwich["lines"][freeLineId] = (
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["points"]][1],
                2,
            )
            freeLineId += 1

            lineRightBack = freeLineId

            sandwich["lines"][freeLineId] = (
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["points"]][2],
                6,
            )
            freeLineId += 1

            lineLeftBack = freeLineId
            sandwich["lines"][freeLineId] = (
                5,
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["points"]][3],
            )
            freeLineId += 1

            replaceLine(
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["lines"]][0],
                [-lineLeftFront, 1, -lineRightFront],
                sandwich,
            )
            replaceLine(
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["lines"]][2],
                [lineRightBack, -5, lineLeftBack],
                sandwich,
            )

            # Remove lineLoop
            sandwich["lineLoops"].pop(
                maxGlueId + layers1.merge[5]["lineLoop"] + maxLayer2Id
            )

            # Remove surface
            sandwich["surfaces"].pop(
                maxGlueId + layers1.merge[5]["surface"] + maxLayer2Id
            )
            sandwich["surfaceLoops"][
                layers1.merge[5]["surfaceLoop"] + maxGlueId + maxLayer2Id
            ].remove(maxGlueId + layers1.merge[5]["surface"] + maxLayer2Id)

            # make lineLoops
            lineLoopLeft = glue.merge[6]["mainLoopGlue"][4] + [
                lineLeftBack,
                -[Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["lines"]][-1],
                -lineLeftFront,
            ]

            lineLoopRight = glue.merge[6]["mainLoopGlue"][3] + [
                -lineRightFront,
                [Id + maxGlueId + maxLayer2Id for Id in layers1.merge[5]["lines"]][1],
                lineRightBack,
            ]

            sandwich["lineLoops"][freeLineLoopId] = lineLoopLeft
            freeLineLoopId += 1

            sandwich["lineLoops"][freeLineLoopId] = lineLoopRight
            freeLineLoopId += 1

            # make surfaces
            surfaceLeft = freeSurfaceId
            sandwich["surfaces"][surfaceLeft] = [freeLineLoopId - 2]
            freeSurfaceId += 1

            surfaceRight = freeSurfaceId
            sandwich["surfaces"][surfaceRight] = [freeLineLoopId - 1]
            freeSurfaceId += 1

            subSurfacesDown = list(
                range(freeSurfaceId, freeSurfaceId + len(glue.merge[6]["subLoops"]))
            )
            for sl in glue.merge[6]["subLoops"]:
                sandwich["surfaces"][freeSurfaceId] = [sl]
                freeSurfaceId += 1

            # make SurfaceLoop
            sandwich["surfaceLoops"][
                layers1.merge[5]["surfaceLoop"] + maxGlueId + maxLayer2Id
            ] = (
                sandwich["surfaceLoops"][
                    layers1.merge[5]["surfaceLoop"] + maxGlueId + maxLayer2Id
                ]
                + [surfaceLeft, surfaceRight]
                + subSurfacesDown
                + [glue.merge[6]["surfaceId"]]
            )

        self.geometry = sandwich

    def shift(self, shiftFactor):
        shiftedPoints = {
            Id + shiftFactor: self.geometry["points"][Id]
            for Id in self.geometry["points"]
        }
        shiftedLines = {
            Id
            + shiftFactor: tuple(
                pId + shiftFactor for pId in self.geometry["lines"][Id]
            )
            for Id in self.geometry["lines"]
        }
        shiftedLineLoops = {
            Id
            + shiftFactor: [
                lId + shiftFactor if lId > 0 else lId - shiftFactor
                for lId in self.geometry["lineLoops"][Id]
            ]
            for Id in self.geometry["lineLoops"]
        }
        shiftedSurfaces = {
            Id
            + shiftFactor: [
                llId + shiftFactor for llId in self.geometry["surfaces"][Id]
            ]
            for Id in self.geometry["surfaces"]
        }
        shiftedSurfaceLoops = {
            Id
            + shiftFactor: [
                sId + shiftFactor for sId in self.geometry["surfaceLoops"][Id]
            ]
            for Id in self.geometry["surfaceLoops"]
        }
        shiftedVolumes = {
            Id
            + shiftFactor: [ssId + shiftFactor for ssId in self.geometry["volumes"][Id]]
            for Id in self.geometry["volumes"]
        }

        self.geometry = {
            "points": shiftedPoints,
            "lines": shiftedLines,
            "lineLoops": shiftedLineLoops,
            "surfaces": shiftedSurfaces,
            "surfaceLoops": shiftedSurfaceLoops,
            "volumes": shiftedVolumes,
        }
        self.merge[33]["points"]["left"] = [
            Id + shiftFactor for Id in self.merge[33]["points"]["left"]
        ]
        self.merge[333]["points"]["left"] = [
            Id + shiftFactor for Id in self.merge[333]["points"]["left"]
        ]
        self.merge[44]["points"]["left"] = [
            Id + shiftFactor for Id in self.merge[44]["points"]["left"]
        ]
        self.merge[444]["points"]["left"] = [
            Id + shiftFactor for Id in self.merge[444]["points"]["left"]
        ]

        self.merge[33]["points"]["right"] = [
            Id + shiftFactor for Id in self.merge[33]["points"]["right"]
        ]
        self.merge[333]["points"]["right"] = [
            Id + shiftFactor for Id in self.merge[333]["points"]["right"]
        ]
        self.merge[44]["points"]["right"] = [
            Id + shiftFactor for Id in self.merge[44]["points"]["right"]
        ]
        self.merge[444]["points"]["right"] = [
            Id + shiftFactor for Id in self.merge[444]["points"]["right"]
        ]

        self.merge[33]["lines"]["right"] = [
            Id + shiftFactor for Id in self.merge[33]["lines"]["right"]
        ]
        self.merge[333]["lines"]["right"] = [
            Id + shiftFactor for Id in self.merge[333]["lines"]["right"]
        ]
        self.merge[44]["lines"]["right"] = [
            Id + shiftFactor for Id in self.merge[44]["lines"]["right"]
        ]
        self.merge[444]["lines"]["right"] = [
            Id + shiftFactor for Id in self.merge[444]["lines"]["right"]
        ]

        self.merge[33]["lines"]["left"] = [
            Id + shiftFactor for Id in self.merge[33]["lines"]["left"]
        ]
        self.merge[333]["lines"]["left"] = [
            Id + shiftFactor for Id in self.merge[333]["lines"]["left"]
        ]
        self.merge[44]["lines"]["left"] = [
            Id + shiftFactor for Id in self.merge[44]["lines"]["left"]
        ]
        self.merge[444]["lines"]["left"] = [
            Id + shiftFactor for Id in self.merge[444]["lines"]["left"]
        ]

        self.merge[33]["lines"]["middle"] = [
            Id + shiftFactor for Id in self.merge[33]["lines"]["middle"]
        ]
        self.merge[333]["lines"]["middle"] = [
            Id + shiftFactor for Id in self.merge[333]["lines"]["middle"]
        ]
        self.merge[44]["lines"]["middle"] = [
            Id + shiftFactor for Id in self.merge[44]["lines"]["middle"]
        ]
        self.merge[444]["lines"]["middle"] = [
            Id + shiftFactor for Id in self.merge[444]["lines"]["middle"]
        ]

        self.merge[33]["lineLoops"] = [
            Id + shiftFactor for Id in self.merge[33]["lineLoops"]
        ]
        self.merge[333]["lineLoops"] = [
            Id + shiftFactor for Id in self.merge[333]["lineLoops"]
        ]
        self.merge[44]["lineLoops"] = [
            Id + shiftFactor for Id in self.merge[44]["lineLoops"]
        ]
        self.merge[444]["lineLoops"] = [
            Id + shiftFactor for Id in self.merge[444]["lineLoops"]
        ]

        self.merge[33]["surfaces"] = [
            Id + shiftFactor for Id in self.merge[33]["surfaces"]
        ]
        self.merge[333]["surfaces"] = [
            Id + shiftFactor for Id in self.merge[333]["surfaces"]
        ]
        self.merge[44]["surfaces"] = [
            Id + shiftFactor for Id in self.merge[44]["surfaces"]
        ]
        self.merge[444]["surfaces"] = [
            Id + shiftFactor for Id in self.merge[444]["surfaces"]
        ]

    def maxId(self):
        return max(
            max(self.geometry["points"].keys()),
            max(self.geometry["lines"].keys()),
            max(self.geometry["lineLoops"].keys()),
        )


class Blade:
    def __init__(self, s1, s2, s3, s4):
        """
            1
        4       2
            3

        Connect layer2 with layer1, clockwise



        """

        self.sandwichUp = s1
        self.sandwichRight = s2
        self.sandwichDown = s3
        self.sandwichLeft = s4
        self.nLayers = s1.nLayers
        lGeometries = [s1.geometry, s2.geometry, s3.geometry, s4.geometry]
        self.geometry = mergeGeometries(lGeometries)

        # ========================================
        # Now make the connections
        # ========================================
        freeLineId = max(self.geometry["lines"].keys()) + 1
        freeLineLoopId = max(self.geometry["lineLoops"].keys()) + 1
        freeSurfaceId = max(self.geometry["surfaces"].keys()) + 1
        freeSurfaceLoopId = max(self.geometry["surfaceLoops"].keys()) + 1
        freeVolumeId = max(self.geometry["volumes"].keys()) + 1

        nLayers = self.nLayers

        lSandwichUp = [
            self.sandwichUp.merge[333],
            self.sandwichRight.merge[333],
            self.sandwichDown.merge[333],
            self.sandwichLeft.merge[333],
        ]
        lSandwichDown = [
            self.sandwichRight.merge[44],
            self.sandwichDown.merge[44],
            self.sandwichLeft.merge[44],
            self.sandwichUp.merge[44],
        ]
        for j in range(4):
            sandwichUp = lSandwichUp[j]
            sandwichDown = lSandwichDown[j]

            # ========================================
            # MAKE LINES
            # ========================================

            linesFront = range(freeLineId, freeLineId + self.nLayers + 1)
            for i in range(nLayers + 1):
                self.geometry["lines"][freeLineId] = (
                    sandwichUp["points"]["left"][i],
                    sandwichDown["points"]["right"][i],
                )
                freeLineId += 1

            linesBack = range(freeLineId, freeLineId + nLayers + 1)
            for i in range(self.nLayers + 1):
                self.geometry["lines"][freeLineId] = (
                    sandwichUp["points"]["right"][i],
                    sandwichDown["points"]["left"][i],
                )
                freeLineId += 1

            # ========================================
            # MAKE LINELOOPS
            # ========================================
            # FRONT
            lineLoopsFront = range(freeLineLoopId, freeLineLoopId + nLayers)
            for i in range(nLayers):
                self.geometry["lineLoops"][freeLineLoopId] = [
                    linesFront[i],
                    -sandwichDown["lines"]["right"][i],
                    -linesFront[i + 1],
                    -sandwichUp["lines"]["left"][i],
                ]
                freeLineLoopId += 1

            # BACK
            lineLoopsBack = range(freeLineLoopId, freeLineLoopId + nLayers)
            for i in range(nLayers):
                self.geometry["lineLoops"][freeLineLoopId] = [
                    -linesBack[i],
                    -sandwichUp["lines"]["right"][i],
                    linesBack[i + 1],
                    -sandwichDown["lines"]["left"][i],
                ]
                freeLineLoopId += 1

            # CROSS SECTION
            lineLoopsCross = range(freeLineLoopId, freeLineLoopId + nLayers + 1)
            for i in range(nLayers + 1):
                self.geometry["lineLoops"][freeLineLoopId] = [
                    linesFront[i],
                    sandwichDown["lines"]["middle"][i],
                    -linesBack[i],
                    -sandwichUp["lines"]["middle"][i],
                ]
                freeLineLoopId += 1

            # ========================================
            # MAKE SURFACES
            # ========================================
            # FRONT
            surfacesFront = range(freeSurfaceId, freeSurfaceId + nLayers)
            for i in range(nLayers):
                self.geometry["surfaces"][freeSurfaceId] = [lineLoopsFront[i]]
                freeSurfaceId += 1
            # BACK
            surfacesBack = range(freeSurfaceId, freeSurfaceId + nLayers)
            for i in range(nLayers):
                self.geometry["surfaces"][freeSurfaceId] = [lineLoopsBack[i]]
                freeSurfaceId += 1
            # CROSS SECTION
            surfacesCross = range(freeSurfaceId, freeSurfaceId + nLayers + 1)
            for i in range(nLayers + 1):
                self.geometry["surfaces"][freeSurfaceId] = [lineLoopsCross[i]]
                freeSurfaceId += 1

            # ========================================
            # MAKE SURFACE LOOPS
            # ========================================

            surfaceLoops = range(freeSurfaceLoopId, freeSurfaceLoopId + nLayers)
            for i in range(nLayers):
                self.geometry["surfaceLoops"][freeSurfaceLoopId] = [
                    surfacesCross[i],
                    surfacesFront[i],
                    surfacesBack[i],
                    surfacesCross[i + 1],
                    sandwichUp["surfaces"][i],
                    sandwichDown["surfaces"][i],
                ]
                freeSurfaceLoopId += 1

            # ========================================
            # MAKE VOLUMES
            # ========================================
            for i in range(nLayers):
                self.geometry["volumes"][freeVolumeId] = [surfaceLoops[i]]
                freeVolumeId += 1


def makeSomeSandwiches(nLayers=5):
    depth = 10
    width = 2

    hightGlue = 0.2
    hightLayer1 = 1
    hightLayer2 = 1

    totalHight = hightGlue + hightLayer1 + hightLayer2
    nLayer1 = nLayers
    nLayer2 = nLayer1

    x = 0
    y = 0
    z = 0

    lcoordsGlue = [
        (x, hightLayer1, z),
        (x + width, hightLayer1, z),
        (x + width, hightLayer1 + hightGlue, z),
        (x, hightLayer1 + hightGlue, z),
    ]

    lcoordLayer1 = [(x, y, z), (x + width, y, z)]

    lcoordLayer2 = [
        (x, y + hightLayer1 + hightGlue, z),
        (x + width, y + hightLayer1 + hightGlue, z),
    ]

    # =================================
    # FIRST SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)
    # for i in range(10):
    #    glue.sampleBubble()
    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich1 = Sandwich(w1, glue, w2)

    # =================================
    # SECOND SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich2 = Sandwich(w1, glue, w2)
    sandwich2.shift(sandwich1.maxId())

    for pointID, pointDict in sandwich2.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (pointCoords[1] + width + 1, -pointCoords[0] - 1, pointCoords[2])
        sandwich2.geometry["points"][pointID]["coords"] = newCoords

    # =================================
    # THIRD SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich3 = Sandwich(w1, glue, w2)
    sandwich3.shift(sandwich2.maxId())

    for pointID, pointDict in sandwich3.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (
            -pointCoords[0] + width,
            -pointCoords[1] - width - 2,
            pointCoords[2],
        )
        sandwich3.geometry["points"][pointID]["coords"] = newCoords

    # =================================
    # FOURTH SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich4 = Sandwich(w1, glue, w2)
    sandwich4.shift(sandwich3.maxId())

    for pointID, pointDict in sandwich4.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (-pointCoords[1] - 1, pointCoords[0] - width - 1, pointCoords[2])
        sandwich4.geometry["points"][pointID]["coords"] = newCoords

    return sandwich1, sandwich2, sandwich3, sandwich4


def makeSomeSandwiches2(nLayers=5):
    depth = 2
    width = 2

    hightGlue = 0.2
    hightLayer1 = 1
    hightLayer2 = 1

    totalHight = hightGlue + hightLayer1 + hightLayer2
    nLayer1 = 2
    nLayer2 = nLayer1

    x = 0
    y = 0
    z = 0

    lcoordsGlue = [
        (x + 0.1, hightLayer1, z),
        (x + width - 0.1, hightLayer1, z),
        (x + width - 0.1, hightLayer1 + hightGlue, z),
        (x + 0.1, hightLayer1 + hightGlue, z),
    ]

    lcoordLayer1 = [(x, y, z), (x + width, y, z)]

    lcoordLayer2 = [
        (x, y + hightLayer1 + hightGlue, z),
        (x + width, y + hightLayer1 + hightGlue, z),
    ]

    # =================================
    # FIRST SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)
    for i in range(4):
        glue.sampleBubble()
    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich1 = Sandwich(w1, glue, w2)

    # =================================
    # SECOND SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)
    for i in range(4):
        glue.sampleBubble()

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich2 = Sandwich(w1, glue, w2)
    sandwich2.shift(sandwich1.maxId())

    for pointID, pointDict in sandwich2.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (pointCoords[1] + width + 1, -pointCoords[0] - 1, pointCoords[2])
        sandwich2.geometry["points"][pointID]["coords"] = newCoords

    # =================================
    # THIRD SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)
    for i in range(4):
        glue.sampleBubble()

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich3 = Sandwich(w1, glue, w2)
    sandwich3.shift(sandwich2.maxId())

    for pointID, pointDict in sandwich3.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (
            -pointCoords[0] + width,
            -pointCoords[1] - width - 2,
            pointCoords[2],
        )
        sandwich3.geometry["points"][pointID]["coords"] = newCoords

    # =================================
    # FOURTH SANDWICH
    # =================================

    glue = Glue(lcoordsGlue, depth)
    for i in range(4):
        glue.sampleBubble()

    glue._finish_geometry()

    w1 = Waffle(
        lcoords=lcoordLayer1,
        depth=depth,
        nLayers=nLayer1,
        hightLayer=hightLayer1 / nLayer1,
    )

    w2 = Waffle(
        lcoords=lcoordLayer2,
        depth=depth,
        nLayers=nLayer2,
        hightLayer=hightLayer2 / nLayer2,
    )

    sandwich4 = Sandwich(w1, glue, w2)
    sandwich4.shift(sandwich3.maxId())

    for pointID, pointDict in sandwich4.geometry["points"].items():
        pointCoords = pointDict["coords"]
        newCoords = (-pointCoords[1] - 1, pointCoords[0] - width - 1, pointCoords[2])
        sandwich4.geometry["points"][pointID]["coords"] = newCoords

    return sandwich1, sandwich2, sandwich3, sandwich4


def mergeGeometries(lGeometries):
    lpoints = [geo["points"] for geo in lGeometries]
    llines = [geo["lines"] for geo in lGeometries]
    llineLoops = [geo["lineLoops"] for geo in lGeometries]
    lsurfaces = [geo["surfaces"] for geo in lGeometries]
    lsurfaceLoops = [geo["surfaceLoops"] for geo in lGeometries]
    lvolumes = [geo["volumes"] for geo in lGeometries]

    geometry = {
        "points": {k: v for x in lpoints for k, v in x.items()},
        "lines": {k: v for x in llines for k, v in x.items()},
        "lineLoops": {k: v for x in llineLoops for k, v in x.items()},
        "surfaces": {k: v for x in lsurfaces for k, v in x.items()},
        "surfaceLoops": {k: v for x in lsurfaceLoops for k, v in x.items()},
        "volumes": {k: v for x in lvolumes for k, v in x.items()},
    }

    return geometry
