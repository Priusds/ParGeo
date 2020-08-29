from .glue import *
from .layers import *
from ..geo_utils.utils import write_geo

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

    assert IdOld > 0
    oldLine = geometry["lines"][IdOld]

    lnewID = list(newLines[0])

    lnewlines = list(newLines[1])

    # Check if new and old line have same orientation
    startOld = oldLine[0]
    startNew = lnewlines[0][0] if lnewID[0] > 0 else lnewlines[0][1]
    sameOrientation = (startOld == startNew)

    if not sameOrientation:
        lnewID = [-l if l is not None else None for l in reversed(lnewID)]
        lnewlines.reverse()
    

    geometry["lines"].pop(IdOld, None)

    for i,newID in enumerate(lnewID): # add new lines
        if newID is None:
            freeLineID = max(geometry["lines"].keys())+1
            
            lnewID[i] = freeLineID
            
            endPoint = geometry["lines"][lnewID[i-1]][1] if lnewID[i-1]>0 else geometry["lines"][-lnewID[i-1]][0]
            if endPoint != lnewlines[i][0]:
                lnewlines[i] = (endPoint, lnewlines[i][0])
            geometry["lines"][freeLineID] = lnewlines[i]
        else:
            geometry["lines"][abs(newID)] = lnewlines[i]



    for ll, llines in geometry["lineLoops"].items():
        #if ll == 555:
           # print("TEST3", llines, IdOld)
        
        if IdOld in llines:
            idx = llines.index(IdOld)
            geometry["lineLoops"][ll] = llines[0:idx] + lnewID + llines[idx+1:]
            return lnewID
           # print(lnewID)
        elif -IdOld in llines:
            idx = llines.index(-IdOld)
            geometry["lineLoops"][ll] = llines[0:idx] + [-l for l in reversed(lnewID)] + llines[idx+1:]   
            return   [-l for l in reversed(lnewID)]
          #  print(geometry["lineLoops"][ll], "neg", IdOld, ll)       
        else:
            pass
    
    


def replaceLineLoop(oldID, newLineLoop, geometry):
    newID = newLineLoop[0]
    newLines = newLineLoop[1]
    
    assert newID > 0

    geometry["lineLoops"].pop(oldID)
    geometry["lineLoops"][newID] = newLines
    for surfaceID, lineLoops in geometry["surfaces"].items():
        if oldID in lineLoops:
            idx = lineLoops.index(oldID)
            geometry["surfaces"][surfaceID] = lineLoops[0:idx] + [newID] + lineLoops[idx+1:]
        elif -oldID in lineLoops:
            idx = lineLoops.index(-oldID)
            geometry["surfaces"][surfaceID] = lineLoops[0:idx] + [-newID] + lineLoops[idx+1:]
        else:
            pass

def replaceSurface(oldID, newSurface, geometry):

    newID = newSurface[0]
    newLineLoops = newSurface[1]

    assert newID > 0

    geometry["surfaces"].pop(oldID)
    geometry["surfaces"][newID] = newLineLoops
    for sLoopID, surfaces in geometry["surfaceLoops"].items():

        if oldID in surfaces:
            idx = surfaces.index(oldID)
            geometry["surfaceLoops"][sLoopID] = surfaces[0:idx] + [newID] + surfaces[idx+1:]
        elif -oldID in surfaces:
            idx = surfaces.index(-oldID)
            geometry["surfaceLoops"][sLoopID] = surfaces[0:idx] + [-newID] + surfaces[idx+1:]
        else:
            pass


def merge(glue, layers):
    assert glue.ready == True
    # First shift indices 
    maxPointId=max(glue.geometry["points"].keys())
    maxLineId=max(glue.geometry["lines"].keys())
    maxLineLoopId=max(glue.geometry["lineLoops"].keys())
    maxSurfaceId=max(glue.geometry["surfaces"].keys())
    maxSurfaceLoopId=max(glue.geometry["surfaceLoops"].keys())
    maxVolumeId=max(glue.geometry["volumes"].keys())


    newPoints = {maxPointId+Id:layers.geometry["points"][Id] for Id in layers.geometry["points"].keys()}
    newLines = {Id+maxLineId:(layers.geometry["lines"][Id][0]+maxPointId,layers.geometry["lines"][Id][1]+maxPointId) for Id in layers.geometry["lines"].keys()}
    newLineLoops = {Id+maxLineLoopId:[i+maxLineId if i > 0 else i - maxLineId for i in layers.geometry["lineLoops"][Id]] for Id in layers.geometry["lineLoops"].keys()}
    newSurfaces = {Id+maxSurfaceId:[i+maxLineLoopId if i > 0 else i - maxLineLoopId for i in layers.geometry["surfaces"][Id]] for  Id in layers.geometry["surfaces"].keys()}
    newSurfaceLoops = {Id+maxSurfaceLoopId:[i+maxSurfaceId if i > 0 else i - maxSurfaceId for i in layers.geometry["surfaceLoops"][Id]] for Id in layers.geometry["surfaceLoops"].keys()}
    newVolumes = {Id+maxVolumeId:[i+maxSurfaceLoopId if i > 0 else i - maxSurfaceLoopId for i in layers.geometry["volumes"][Id]] for Id in layers.geometry["volumes"].keys()}
    
    shiftyLayers = {
        "points":newPoints,
        "lines":newLines,
        "lineLoops":newLineLoops,
        "surfaces":newSurfaces,
        "surfaceLoops":newSurfaceLoops,
        "volumes":newVolumes
    }
    LayerPart = shiftyLayers
    GluePart = glue.geometry
    TotalGeometry = {
        "points":{**LayerPart["points"], ** GluePart["points"]},
        "lines":{**LayerPart["lines"], ** GluePart["lines"]},
        "lineLoops":{**LayerPart["lineLoops"], ** GluePart["lineLoops"]},
        "surfaces":{**LayerPart["surfaces"], ** GluePart["surfaces"]},
        "surfaceLoops":{**LayerPart["surfaceLoops"], ** GluePart["surfaceLoops"]},
        "volumes":{**LayerPart["volumes"], ** GluePart["volumes"]}
    }

    #freePointId=max(newPoints.keys())+1
    #freeLineId=max(newLines.keys())+1
    freeLineLoopId=max(newLineLoops.keys())+1
    freeSurfaceId=max(newSurfaces.keys())+1
    #freeSurfaceLoopId=max(newSurfaceLoops.keys())+1
    #freeVolumeId=max(newVolumes.keys())+1


    intersectionFacet = 5 # From glue perspective ==>> layers is above
    interSectionFacetLayer = 6
    case = 1
    if case == 1:
        # Glue and Layer have same width everywhere

        # Replace points
        lnewPointID = [4,3,7,8] 
        for point,newID in zip(layers.PointsDown, lnewPointID):
            
            replacePoint(point+maxPointId, (newID, GluePart["points"][newID]), TotalGeometry)
        
        # Replace lines
        print("TEST:", glue.faces[5]["outerLoop"])
        newEdges = []
        for i in range(1,5):
            
            newEdge = replaceLine([l +maxLineId for l in layers.LinesDown][i-1], list(zip(*glue.faces[intersectionFacet]["outerLoop"][i]["edgeLines"])), TotalGeometry)
            newEdges= newEdges + newEdge
        
        # Make new LineLoop
        
        #TotalGeometry["lineLoops"][freeLineLoopId] = newEdges
        
        
        newLoop = (freeLineLoopId, newEdges)
      #  print(layers.LineLoopDown + maxLineLoopId)
      #  print(newLoop)
        replaceLineLoop(layers.LineLoopDown+maxLineLoopId, newLoop,TotalGeometry)
        #print("MainLoop, ",  glue.faces[intersectionFacet]["merge"]["mainLoop"])
        # Make new Surface and replace Old one
        lineLoops = [freeLineLoopId] + glue.faces[intersectionFacet]["innerLoops"]

        newSurface = (freeSurfaceId, lineLoops)
        freeSurfaceId += 1
        #print(layers.SurfaceDown+ maxSurfaceId)
        replaceSurface(layers.SurfaceDown + maxSurfaceId, newSurface, TotalGeometry)

        # Add hole surfaces

        for innerLoop in glue.faces[intersectionFacet]["innerLoops"]:
            TotalGeometry["surfaces"][freeSurfaceId] = [innerLoop]
            freeSurfaceId += 1
        nHoles = len(glue.faces[intersectionFacet]["innerLoops"])
        newSurfaceIDs = list(range(freeSurfaceId-nHoles, freeSurfaceId))
        TotalGeometry["surfaceLoops"][3] = TotalGeometry["surfaceLoops"][3] + newSurfaceIDs
        #print(LayerPart["surfaceLoops"])
        #print([glue.faces[5]["outerLoop"][i]["lines"][0]  for i in glue.faces[5]["outerLoop"].keys() if i is int])
       # mainBubbleLoop = glue.faces[intersectionFacet]["merge"]["mainLoop"]
       # lineLoopsBubblOutList = []
       # lInnerBubbles = glue.faces[intersectionFacet]["innerLoops"]
        
      #  for i in range(1,5):
      #      lineLoopsBubblOutList = lineLoopsBubblOutList + glue.faces[intersectionFacet]["outerLoop"][i]["lineLoops"]
        
        # Iterate over lineLoopsBubbl create new lines ==>> new Loop
       # for loop in lineLoopsBubblOutList:
       #     firstLine = loop[0]
       #     lastLine = loop[-1]
       #     firstPoint = glue.geometry["lines"][firstLine][0] if firstLine > 0 else glue.geometry["lines"][-firstLine][1]
       #     lastPoint = glue.geometry["lines"][lastLine][1] if firstLine > 0 else glue.geometry["lines"][-lastLine][0]
       #     IntersectionPart["lines"][freeLineId] = (lastPoint, firstPoint)
       #     IntersectionPart["lineLoops"][freeLineLoopId] =  loop + [freeLineId]
       #     freeLineId += 1
       #     freeLineLoopId+=1
        
        # remove intersection Surface from layer part ==>> change SurfaceLoop and Volume
        
        #LayerPart["surfaces"].pop(layers.SurfaceDown)



    else: 
        raise NotImplementedError

    

    write_geo(TotalGeometry, "Total")
    #write_geo(GluePart, "glucamole")
    #write_geo(LayerPart, "waffel")