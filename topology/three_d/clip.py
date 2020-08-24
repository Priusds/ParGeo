from abc import ABC, abstractmethod
from .utils import locate, normalize
import math


import mystic as my
# tools
from mystic.monitors import VerboseMonitor
from mystic.math.measures import mean, impose_mean
from mystic.math import almostEqual
from mystic.penalty import quadratic_equality, quadratic_inequality

import numpy as np


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



class NumericPointProjection(PointProjection):
    def __init__(self, projectionPlane, trafo):
        self.normalVector = normalize(projectionPlane[0])
        self.supportVector = projectionPlane[1]
        self.trafo = trafo

    def apply(self, point):

        P = np.array(point)

        f = lambda p : np.linalg.norm(self.trafo(p) - P)**2
        constraint = lambda p : (self.trafo(p) - np.array(self.supportVector)).transpose().dot(np.array(self.normalVector))


        x0 = np.array([0.5*np.pi, np.pi])


        def penalty01(x):  # <= 0.0
            return x[0]
        def penalty02(x):  # <= 0.0
            return x[0] - np.pi
        def penalty11(x):  # <= 0.0
            return x[1]
        def penalty12(x):  # <= 0.0
            return x[1] - 2*np.pi
        
        def penalty2(u):  # ==0
            return constraint(u) # np.abs()**2 - 0.0001

        #@quadratic_inequality(penalty01, k=1e12)
        #@quadratic_inequality(penalty02, k=1e12)
        #@quadratic_inequality(penalty11, k=1e12)
        #@quadratic_inequality(penalty12, k=1e12)
        @quadratic_inequality(penalty2, k=1e12)
        def penalty(x):
            return 0.0

        solver = my.solvers.fmin
        verbose = False
        if verbose:
            mon = my.monitors.VerboseMonitor(100)
        else:
            mon = None

        kwds = dict(disp=verbose, full_output=True, #itermon=mon,
                    xtol=1e-12, ftol=1e-12, maxfun=10000, maxiter=100)
        result = solver(f, x0, penalty=penalty, **kwds)
        print("START WITH P = ", point)
        print("x0 : ", np.array([0.5*np.pi, np.pi]))
        print("u* = {v} -> g(u*) = {w}".format(v=result[0], w=constraint(result[0])))

        res =  tuple(self.trafo(result[0]))
        print(" new point : ", res)
        return res





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


