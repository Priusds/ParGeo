from abc import ABC, abstractmethod
from .utils import locate, normalize
import math

class PointProjection(ABC):
    def __init__(self, projectionPlane):
        """
            defines the hyperplane where we project the data.
        """
        pass
    def apply(self, point):
        """
            point = {"id","coords"}
            self.apply(point) changes point["coords"]

            projects the point on the projection area.
        """
        pass

class EuclideanProjection(PointProjection):
    def __init__(self, projectionPlane):
        self.normalVector = normalize(projectionPlane[0])
        self.supportVector = projectionPlane[1]

    def apply(self, point):
        center_point = [point[i] - self.supportVector[i] for i in range(3)]
        normalPart = sum([center_point[i]*self.normalVector[i] for i in range(3)])
        project_point = [center_point[i] - normalPart*self.normalVector[i] for i in range(3)]
        return tuple([project_point[i] + self.supportVector[i] for i in range(3)])

class ClippingRule(ABC):
    def __init__(self, pointProjection):
        """
            contains the main informations needed to decide wich points should be projected

            changes 
        """
        raise NotImplementedError

    def apply(self, ltriangle, lline, lpoint, geometry):
        raise NotImplementedError

class EasyClip(ClippingRule):
    def __init__(self, PointProjection):
        self.normalVector = PointProjection.normalVector
        self.supportVector = PointProjection.supportVector
        self.projectPoint = PointProjection.apply
    
    def apply(self, ltriangle, lline, lpoint, geometry):
        """
            project only points from outside on boundary, e.g. only change some coords of lpoint
        """
        points = dict()
        for Id in lpoint:
            coords = geometry["points"][Id]["coords"]
            if not locate(coords, self.normalVector, self.supportVector):
                # point is outside
                coords = self.projectPoint(coords)
            points[Id] = {"coords":coords, "lc":geometry["points"][Id]["lc"]}
        return ltriangle, lline, points
