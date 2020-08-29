from ..geo_utils.utils import write_geo


class Waffle:

    def __init__(self, lcoords, depth, nLayers,hightLayer):
        # lcoords = [c1, c2]
        self.geometry = {
            "points":dict(),
            "lines":dict(),
            "lineLoops":dict(),
            "surfaces":dict(),
            "surfaceLoops":{}, # 1 is outside, 0 is inside
            "volumes":{}
        }
        self.min_y = min(lcoords[0][1], lcoords[1][1])
        self.layerHights = [n*hightLayer for n in range(nLayers+1)]
        self.max_y = self.layerHights[-1]+hightLayer

        self.max_x = lcoords[1][0]
        self.max_z = lcoords[1][2]
        self.min_x = lcoords[0][0]
        
        self.min_z = -depth

        

        self.freePointId = 1
        self.freeLineId = 1
        self.freeLineLoopId = 1
        self.freeSurfaceId = 1
        self.freeSurfaceLoopId = 1
        self.freeVolumeId = 1

        # =====     POINTS      =======

        # Front:
        # x and z coords are constant
        firstPointFrontLeftId = self.freePointId
        firstPointFrontRightId = firstPointFrontLeftId + nLayers+1
        xLeft = self.min_x
        xRight = self.max_x
        z = self.max_z 

        for i in range(nLayers+1):
            yLeft =lcoords[0][1]+ self.layerHights[i]
            yRight = lcoords[1][1]+self.layerHights[i]
            Id = self.freePointId
            pointLeft = {"coords":(xLeft, yLeft, z), "lc":None}
            pointRight = {"coords":(xRight, yRight, z), "lc":None}
            self.freePointId += 1
            self.geometry["points"][Id]=pointLeft
            self.geometry["points"][Id+nLayers+1]=pointRight

        self.freePointId+=nLayers+1

        # Back:
        firstPointBackLeftId = self.freePointId
        firstPointBackRightId = firstPointBackLeftId + nLayers+1
        xLeft = self.max_x
        xRight = self.min_x
        z = self.min_z 

        for i in range(nLayers+1):
            yLeft =lcoords[0][1]+ self.layerHights[i]
            yRight = lcoords[1][1]+self.layerHights[i]
            Id = self.freePointId
            pointLeft = {"coords":(xLeft, yRight, z), "lc":None}
            pointRight = {"coords":(xRight, yLeft, z), "lc":None}
            self.freePointId += 1
            self.geometry["points"][Id]=pointLeft
            self.geometry["points"][Id+nLayers+1]=pointRight
        self.freePointId+=nLayers+1


        # ======    LINES     ========

        # Front
        firstLineFrontLeftVId = self.freeLineId
        firstLineFrontRightVId = firstLineFrontLeftVId + nLayers
        for i in range(nLayers):
            IdLeft = self.freeLineId
            startLeft = firstPointFrontLeftId + i
            endLeft = firstPointFrontLeftId + i +1
            self.geometry["lines"][IdLeft] = (endLeft, startLeft)

            IdRight = IdLeft + nLayers
            startRight = startLeft + nLayers+1
            endRight = startRight +1
            self.geometry["lines"][IdRight] = (startRight,endRight)
            self.freeLineId+=1

        self.freeLineId += nLayers
        firstLineFrontHId = self.freeLineId
        for i in range(nLayers+1):
            start = firstPointFrontLeftId+i
            end = start + nLayers+1
            self.geometry["lines"][self.freeLineId]=(start, end)
            self.freeLineId+=1

        # Back
        firstLineBackLeftVId = self.freeLineId
        firstLineBackRightVId = firstLineBackLeftVId + nLayers
        for i in range(nLayers):
            IdLeft = self.freeLineId
            startLeft = firstPointBackLeftId + i
            endLeft = firstPointBackLeftId + i +1
            self.geometry["lines"][IdLeft] = (endLeft, startLeft)

            IdRight = IdLeft + nLayers
            startRight = startLeft + nLayers+1
            endRight = startRight +1
            self.geometry["lines"][IdRight] = (startRight,endRight)
            self.freeLineId+=1
            
        self.freeLineId += nLayers
        firstLineBackHId = self.freeLineId
        for i in range(nLayers+1):
            start = firstPointBackLeftId+i
            end = start + nLayers+1
            self.geometry["lines"][self.freeLineId]=(start, end)
            self.freeLineId+=1
        
        # Connection Front Back
        firstConnectionLeft = self.freeLineId
        for i in range(nLayers+1):
            start = firstPointFrontLeftId+i
            end = firstPointBackRightId+i
            self.geometry["lines"][self.freeLineId] = (start, end)
            self.freeLineId+=1
        firstConnectionRight = self.freeLineId
        for i in range(nLayers+1):
            start = firstPointFrontRightId+i
            end = firstPointBackLeftId+i
            self.geometry["lines"][self.freeLineId] = (start, end)
            self.freeLineId+=1

        # ===== LINELOOPS   =======
        # LineLoopsFront
        firstLineLoopFront = self.freeLineLoopId
        for i in range(nLayers):
            l1 = firstLineFrontHId+i
            l3 = -(firstLineFrontHId+1+i)
            l2 = firstLineFrontRightVId+i
            l4 = firstLineFrontLeftVId+i
            self.geometry["lineLoops"][self.freeLineLoopId]=[l1,l2,l3,l4]
            self.freeLineLoopId+=1

        # LineLoopsBack
        firstLineLoopBack = self.freeLineLoopId
        for i in range(nLayers):
            l1 = firstLineBackHId+i
            l3 = -(firstLineBackHId+1+i)
            l2 = firstLineBackRightVId+i
            l4 = firstLineBackLeftVId+i
            self.geometry["lineLoops"][self.freeLineLoopId]=[l1,l2,l3,l4]
            self.freeLineLoopId+=1

        # LineLoopsLeft
        firstLineLoopLeft = self.freeLineLoopId
        for i in range(nLayers):
            l1 = -(firstConnectionLeft+i)
            l3 = -l1 +1
            l2 = -(firstLineFrontLeftVId+i)
            l4 = -(firstLineBackRightVId+i)
            self.geometry["lineLoops"][self.freeLineLoopId]=[l1,l2,l3,l4]
            self.freeLineLoopId+=1

        # LineLoopsRight
        firstLineLoopRight = self.freeLineLoopId
        for i in range(nLayers):
            l1 = firstConnectionRight+i
            l3 = -(l1 +1)
            l2 = -(firstLineFrontRightVId+i)
            l4 = -(firstLineBackLeftVId+i)
            self.geometry["lineLoops"][self.freeLineLoopId]=[l1,l2,l3,l4]
            self.freeLineLoopId+=1


        # LineLoopCrossSection
        firstLineLoopCrossSection = self.freeLineLoopId
        for i in range(nLayers+1):
            l1 = firstLineFrontHId+i
            l2 = firstConnectionRight+i
            l4 = -(firstConnectionLeft+i)
            l3 = firstLineBackHId+i

            self.geometry["lineLoops"][self.freeLineLoopId]=[l1,l2,l3,l4]
            self.freeLineLoopId+=1
        # ====== SURFACES ==========
        # Front Facet
        for i in range(nLayers):
            self.geometry["surfaces"][self.freeSurfaceId]=[firstLineLoopFront+i]
            self.freeSurfaceId+=1
        # Back Facet
        for i in range(nLayers):
            self.geometry["surfaces"][self.freeSurfaceId]=[firstLineLoopBack+i]
            self.freeSurfaceId+=1
        # Left Facet
        for i in range(nLayers):
            self.geometry["surfaces"][self.freeSurfaceId]=[firstLineLoopLeft+i]
            self.freeSurfaceId+=1
        # right Facet
        for i in range(nLayers):
            self.geometry["surfaces"][self.freeSurfaceId]=[firstLineLoopRight+i]
            self.freeSurfaceId+=1
        
        # Cross Sections from bottom to top
        firstCrossSurface = self.freeSurfaceId
        for i in range(nLayers+1):
            self.geometry["surfaces"][self.freeSurfaceId]=[firstLineLoopCrossSection+i]
            self.freeSurfaceId+=1

        # ======= SURFACELOOPS ========
        firstSurfaceLoop = self.freeSurfaceLoopId
        for i in range(nLayers):
            fdown = firstLineLoopCrossSection+i
            fup = firstLineLoopCrossSection+i+1
            fright = firstLineLoopRight+i
            fleft = firstLineLoopLeft+i
            ffront = firstLineLoopFront+i
            fback = firstLineLoopBack+i
            surfaceLoop = [ffront, fback, fright, fleft, fup, fdown]
            self.geometry["surfaceLoops"][self.freeSurfaceLoopId] = surfaceLoop
            self.freeSurfaceLoopId+=1
        
        # ========== VOLUMES ==========0
        for i in range(nLayers):
            self.geometry["volumes"][self.freeVolumeId]=[firstSurfaceLoop+i]
            self.freeVolumeId+=1

        # For merge with Bubbles
        self.PointsDown = [firstPointFrontLeftId, firstPointFrontRightId, firstPointBackLeftId, firstPointBackRightId]
        self.LinesDown = [firstLineFrontHId, firstConnectionRight,firstLineBackHId,firstConnectionLeft]
        self.LineLoopDown = firstLineLoopCrossSection
        self.LineLoopUp = self.LineLoopDown + nLayers + 1
        self.SurfaceDown = firstCrossSurface
        self.SurfaceUp = firstCrossSurface + nLayers + 1
        self.SurfaceLoopDown = firstSurfaceLoop
        self.SurfaceLoopUp = firstSurfaceLoop + nLayers
        #VolumeDown = self.freeVolumeId - nLayers
        #VolumeUp = self.freeVolumeId-1

    def to_geo(self, file_name):
        write_geo(self.geometry, file_name)