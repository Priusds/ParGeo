from abc import ABC, abstractmethod
from .star_geometry import Ellipse


class Logic(ABC):
    """
    Given some subdomain, (rect, hole) where the rect is described by (ax, ay), 
    (since hole and rect have midpoint (0,0) and the rect is not rotated). 
    This class (the subclasses) describes the logic which reference mesh
    to choose. """
    def __init__(self):
        pass
    @abstractmethod
    def number_boundary_points(self, ax,ay):
        raise NotImplementedError

    @abstractmethod 
    def get_reference(self, ax, ay, hole):
        raise NotImplementedError

class Bloed(Logic):
    def __init__(self, c = 3, holerefs = 80):
        self.c = c
        self.refs = holerefs

    def number_boundary_points(self, ax, ay, macro):
        #returns nx, ny
        return self.c, self.c
        
    def get_reference(self,ax, ay, hole):
        refs = self.refs
        r = .95
        ax_ref = 1
        ay_ref = 1
        return [Ellipse((0,0), (r,r)), refs, ax_ref, ay_ref]

    
class CoarseFineLogic(Logic):
    def __init__(self, ccoarse = 1, cfine = 3, holerefs = 80):
        self.c = ccoarse

        self.cfine = cfine
        self.refs = holerefs

    def number_boundary_points(self, ax, ay, macro):
        #returns nx, ny
        if macro == True:
            return self.c, self.c
        else:
            return self.cfine, self.cfine
        
    def get_reference(self,ax, ay, hole):
        refs = self.refs
        r = .95
        ax_ref = 1
        ay_ref = 1
        return [Ellipse((0,0), (r,r)), refs, ax_ref, ay_ref]


class DebugLogic(Logic):
    def __init__(self, c = 3):
        self.c = c

    def number_boundary_points(self, ax, ay):
        #returns nx, ny
        return self.c, self.c
        
    def get_reference(self,ax, ay, hole):
        refs = 4
        r = .5#5
        ax_ref = 1
        ay_ref = 1
        return [Ellipse((0,0), (r,r)), refs, ax_ref, ay_ref]
