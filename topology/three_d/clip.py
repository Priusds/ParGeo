from abc import ABC, abstractmethod


class PointProjection(ABC):
    def __init__(self, projectionPlane):
        """
            defines the hyperplane where we project the data.
        """
        pass
    def apply(point):
        """
            point = {"id","coords"}
            self.apply(point) changes point["coords"]

            projects the point on the projection area.
        """
        pass


class ClippingRule(ABC):
    def __init__(self, pointProjection, triangleLoop, vertices, geo):
        """
            contains the main informations needed to decide wich points should be projected

            changes 
        """
        raise NotImplementedError