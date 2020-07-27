from .utils import get_rect, angle, rotate_clockwise, rotate_counterclockwise, rot_matrix, correct
from abc import ABC, abstractmethod
import numpy as np
import math
from scipy.optimize import minimize as scipy_minimize
import matplotlib.pyplot as plt
cos = math.cos
sin = math.sin
from shapely.geometry import LineString



class StarGeometry(ABC):
    
    def __init__(self, type_, midpoint):
        self.type = type_
        self.midpoint = midpoint

    @abstractmethod
    def discretize_hole(self):
        """ returns a list of points (vertices) in ccw order"""
        raise NotImplementedError

    @abstractmethod
    def to_dict(self):
        """returns the hole informations in dict format"""
        raise NotImplementedError
    
    @abstractmethod
    def get_radius(self, phi):
        """Returns the "radius" in direction phi."""
        raise NotImplementedError

    @abstractmethod
    def rotate(self, phi):
        """Returns a new geometry rotated with angle phi"""
        raise NotImplementedError

    def on_boundary(self, point, tol = 1e-5):
        """Returns True if a point is on the boundary otherwise False."""

        point_centerd = [point[i] - self.midpoint[i] for i in range(2)]
        distance = sum([x**2 for x in point_centerd])
        radius = self.get_radius(angle(point_centerd))

        return abs(radius**2 - distance) < tol

    def info(self):
        """ prints informations about the hole"""
        print(self.to_dict())

    def plot(self, nrefs = 200):
        data = self.discretize_hole(nrefs)
        x,y = list(zip(*data))
        x = list(x)
        y = list(y)
        x.append(x[0])
        y.append(y[0])
        plt.plot(x,y, '-')
    

    @classmethod
    def from_dict(cls, dict):
        if dict['type'] not in LOOKUP_HOLES:
            raise NotImplementedError("from_dict method is not implemented in subclass {0}".format(dict['type']))
        return LOOKUP_HOLES[dict['type']].from_dict(dict)


#=================================SUBCLASSES====================================================



#==========================================================
#=======================____CIRCLE CLASS______=============
#==========================================================
class Circle(StarGeometry):

    def __init__(self, midpoint, radius):
        super().__init__("circle", midpoint)
        self.radius = radius

    @classmethod
    def from_dict(cls, dict):
        return Circle(dict["midpoint"], dict["radius"])

    def to_dict(self):
        return {'type':self.type, 'midpoint': self.midpoint, 'radius':self.radius}

    def discretize_hole(self, refs, coordinates = "euclidean", rwidth = 0., ralphashift = 0.):
        if coordinates == "euclidean":
            return discretize_ellipse([self.midpoint, (self.radius * ( 1 + rwidth), self.radius*(1+rwidth)), 0], refs, ralphashift)
        elif coordinates == "polar":
            return discretize_ellipse_polar([self.midpoint, (self.radius, self.radius), 0], refs)
        else:
            print("error")
            exit()


    def get_radius(self,angle):
        # This method is used for the transform mesh method
        return self.radius

    def rotate(self, angle):
        return self


#==========================================================
#======================____ELLIPSE CLASS______=============
#==========================================================
class Ellipse(StarGeometry):
    def __init__(self, midpoint, axis, angle = 0):
        super().__init__("ellipse", midpoint)
        self.axis = axis
        self.angle = angle

    def to_dict(self):
        return {'type':self.type, 'midpoint':self.midpoint, 'axis': self.axis, 'angle':self.angle}
    
    @classmethod
    def from_dict(cls, dict):
        return Ellipse(dict['midpoint'], dict['axis'], dict['angle'])

    def discretize_hole(self, refs, coordinates = "euclidean", rwidth = 0., ralphashift = 0.):
        if coordinates == "euclidean":
            return discretize_ellipse([self.midpoint, (self.axis[0]*(1+rwidth), self.axis[1]*(1+rwidth)), self.angle], refs, ralphashift =ralphashift)
        elif coordinates == "polar":
            return discretize_ellipse_polar([self.midpoint, self.axis, self.angle], refs)
        else:
            print("error")
            exit()

    def get_radius(self, theta):
        # This method is used for the transform mesh method
        ax, ay = self.axis
        return 1/np.sqrt(np.cos(theta-self.angle)**2 / ax**2 + np.sin(theta - self.angle)**2/ay**2)

    def rotate(self, angle):
        new_angle = correct(self.angle + angle)
        
        return Ellipse(self.midpoint, self.axis, new_angle)


#==========================================================
#======================____STELLAR CLASS______=============
#==========================================================
class Stellar(StarGeometry):
    """
    coeffient must be of type numpy.ndarray
    For smooth forms:
    coeffients = np.random.uniform(-k, k, (j, 2))
    where k in [.1, .2] and j = 1 / k

     """
    def __init__(self, midpoint, radius, coefficient = None):
        super().__init__("stellar", midpoint)
        self.coefficient = np.random.uniform(-.2, .2, (5, 2)) if coefficient is None else coefficient
        self.radius = radius

    @classmethod
    def from_dict(cls, dict):
        if 'coefficient' not in dict:
            raise ValueError("Not possible to initialize stellar hole with using from_dict without coefficients.")
        return Stellar(dict['midpoint'], dict['radius'], np.array(dict['coefficient']))

    def to_dict(self):
        return {'type': self.type, 'midpoint': self.midpoint, 'coefficient':self.coefficient.tolist(), 'radius': self.radius}

    def discretize_hole(self, refs, coordinates = "euclidean", rwidth = 0., ralphashift = 0.):
        return discetize_stellar_polygon(self.midpoint, self.radius* (1+rwidth) , self.coefficient, refs, coordinates, ralphashift )

    
    def get_radius(self, angle):
        # This method is used for the transform mesh method
        return trigonometric_function(self.coefficient, self.radius)(angle)

    def rotate(self, angle ,refs = 80):
        # no analitycal solution 
        new_angle = correct(angle)
        print("angle: ", new_angle)
        return fit_stellar(self, self.midpoint, refs, new_angle)
    
    def flip(self, ncoeff):
        # rotates of 90 degrees

        return flip_fit(self, ncoeff=ncoeff)


    def inflate(self, eps):
        return Stellar(self.midpoint, self.radius*(1+eps), self.coefficient)

    def check_valid(self, npoints = 30):
        return check_valid_midpoint(self,npoints)
#==========================================================
#===================____RECTANGULAR CLASS______============
#==========================================================
class Rectangular(StarGeometry):
    # TODO: add constraint that width > height and 0 <= angle < pi. Maybe?!
    def __init__(self, angle, width, height, midpoint):
        super().__init__("rectangular", midpoint)
        self.angle = angle
        self.width = width
        self.height = height
        self.midpoint = midpoint
        self.boundary = get_rect((angle, 0.5*width, 0.5*height, midpoint))
    
    @classmethod
    def from_dict(cls, dict):
        return Rectangular(dict['angle'], dict['width'], dict['height'], dict['midpoint'])

    def to_dict(self):
        return {'type': self.type, 'midpoint': self.midpoint, 'angle': self.angle, 'width': self.width, 'height': self.height, 'boundary': self.boundary}

    def discretize_hole(self, refs):
        return self.boundary

    def get_radius(self, phi):
        """
        Solve equation: ||c*(cos(phi), sin(phi))||= 1 after c and return cs. 
        The rectangle norm is a weighted infinity norm. 
        ||p|| = max(|p[0]|/ax , |p[1]|/ay).
        """
        p = rotate_clockwise((cos(phi), sin(phi)), self.angle)
        return 1/max(abs(p[0])/(self.width/2), abs(p[1])/(self.height/2))
    
    def rotate(self, angle):
        new_angle = correct(self.angle + angle)
        return Rectangular(new_angle, self.width, self.height, self.midpoint)


#==========================================================
#=================____POLYRECTANGULAR CLASS______==========
#==========================================================
class PolyRectangular(Rectangular):
    """
    The only difference to Rectangular is the boundary. 
    On the edge parallel to x-axis are nx boundary inner points and on y-axis edge there are ny.
    """
    def __init__(self, angle, width, height, midpoint, specific_boundary):
        super().__init__(angle, width, height, midpoint)
        self.type = "poly_rectangular"
        self.boundary = specific_boundary


    @classmethod
    def specify_boundary_nx_ny(cls, rect, nx, ny):
        boundary_vertices = []

        for i in range(4):
            p_i = np.array(rect.boundary[i])
            p_next = np.array(rect.boundary[(i+1) % 4])
            n = (nx,ny)[i%2]

            for j in range(n+1):
                k = (j)/(n+1)
                q = (1-k)*p_i + k*p_next
                boundary_vertices.append(q.tolist())

        return PolyRectangular(rect.angle, rect.width, rect.height, rect.midpoint, boundary_vertices)

    @classmethod
    def from_dict(cls, dict):
        return PolyRectangular(dict['angle'], dict['width'], dict['height'], dict['midpoint'], dict['boundary'])

    def rotate(self, angle):
        new_angle = correct(self.angle + angle)
        return PolyRectangular(new_angle, self.width, self.height, self.midpoint, self.boundary)

LOOKUP_HOLES = {"circle":Circle,
                "stellar":Stellar,
                "ellipse":Ellipse,
                "rectangular":Rectangular,
                "poly_rectangular": PolyRectangular}


def discretize_ellipse_polar(ellipse, refs):
    m, rads, alpha = ellipse
    a,b = rads

    if a >= b:
        angles = [(i / refs) * 2 * np.pi for i in range(refs)]
        eps = math.sqrt(a**2 - b**2)/a
        rad = lambda x: 1/math.sqrt(cos(x - alpha)**2/a**2 + sin(x - alpha)**2/b**2)
        radii = [rad(w) for w in angles]
        return list(zip(radii, angles)), m
    else:
        alpha = alpha + np.pi/2
        c = a
        a = b
        b = c
        angles = [(i / refs) * 2 * np.pi for i in range(refs)]
        eps = math.sqrt(a**2 - b**2)/a
        rad = lambda x: 1/math.sqrt(cos(x - alpha)**2/a**2 + sin(x - alpha)**2/b**2)
        radii = [rad(w) for w in angles]
        return list(zip(radii, angles)), m

def discretize_ellipse(ellipse, refs, ralphashift):

    """
    discretize the ellipse
    :param ellipse: [m, rads, alpha], n
    :return: list of points in ccw
    """

    print("start with ellipse with ", ralphashift)
    print("-----")
    m, rads, alpha = ellipse
    X = np.linspace(0, 2*np.pi, refs, endpoint = False)
    angles = (1-ralphashift)*X + ralphashift*np.concatenate([X[1:], [2*np.pi]])

    circle_points = np.array([
        np.cos(angles),
        np.sin(angles)
    ])

    A = np.array([[rads[0], 0], [0, rads[1]]])
    R = rot_matrix(alpha)
    return ((R @ A @ circle_points).transpose() + m).tolist()


def discetize_stellar_polygon(midpoint, radius, coefficients, refs, coordinates = "euclidean", ralphashift = 0):
    X = np.linspace(0, 2*np.pi, refs, endpoint = False)
    X = (1-ralphashift)*X + ralphashift*np.concatenate([X[1:], [2*np.pi]])
    f = trigonometric_function(coefficients, radius)
    rad = [f(x) for x in X]
    if coordinates == "euclidean":
        return [(r*np.cos(omega) + midpoint[0], r*np.sin(omega) + midpoint[1]) for r, omega in zip(rad, X)]
    else:
        return list(zip(rad,X)),midpoint

def trigonometric_function(coefficients, radius):
    # basis: a_k*cos(k*x) + b_k*sin(k*x)
    """
    Returns a trigonometric functions: angle -> radius

    this method is used to describe the boundary of the stellar hole

    coefficients should be numpy array of shape n x 2"""

    if not type(coefficients) is  np.ndarray:
        raise ValueError("coefficients should be a numpy array")
    if len(coefficients.shape) != 2 or coefficients.shape[1] != 2:
        raise ValueError("coefficients should be a numpy array of shape n x 2")

    n = len(coefficients)
    #return lambda x: 1+ radius * sum([(coefficients[k-1][0]*np.cos(2*k*x) + coefficients[k-1][1]*np.sin(2*k*x)) for k in range(1, n+1)])
    return lambda x: radius * np.exp( sum([1./(k**0.8+1) * (coefficients[k-1][0]*np.cos(k*x) + coefficients[k-1][1]*np.sin(k*x)) for k in range(1, n+1)]))
    #return lambda x: radius * np.exp(sum(coefficients[0])/2 + sum([coefficients[k][0]*np.cos(k*x) + coefficients[k][1]*np.sin(k*x) for k in range(1, n)]))
    # TODO replace np.exp(sum(coefficients[0])/2 with some value
    #return lambda x: radius * np.exp(sum(coefficients[0])/2 + sum([coefficients[k][0]*np.cos(k*x) + coefficients[k][1]*np.sin(k*x) for k in range(1, n)]))


def fit(points, midpoint, parametrisation, initial_value):
    points = np.array(points)
    centered_points = points - midpoint
    n_points = points.shape[0]
    angles = [math.atan2(p[1],p[0]) for p in centered_points]

    p_y = np.array([np.sqrt(centered_points[i][0] ** 2 + centered_points[i][1] ** 2) for i in range(n_points)])
    

    def fun1(x):
        f_y =  np.array([parametrisation(x)(phi) for phi in angles])  #x[0] * np.exp(cos_k @ x[1:m] + sin_k @ x[m:])
        return np.sum((f_y - p_y) ** 2)

    res = scipy_minimize(fun1, initial_value, method='SLSQP')
    print(res.mess)
    return res.x


def flip_fit(hole, refs = 80, ncoeff = 20):
    midpoint = hole.midpoint
    hole.midpoint = (0,0)
    points = hole.discretize_hole(refs)
    R = rot_matrix(np.pi/2)
    rotated_points = np.array([np.array(p) for p in points]) @ R.transpose()
   
    n_points = len(points)
    angles = [math.atan2(p[1],p[0]) for p in rotated_points]
    

    m = int(ncoeff/ 2)

    cos_k = np.array([[np.cos(j * phi)/(j**0.8 + 1) for j in range(1, m+1)] for phi in angles])
    sin_k = np.array([[np.sin(j * phi)/(j**0.8 + 1) for j in range(1, m+1)] for phi in angles])

    p_y = np.array([np.sqrt(rotated_points[i][0] ** 2 + rotated_points[i][1] ** 2) for i in range(n_points)])

   # print(jac(np.zeros(N + 1), p_y, cos_k, sin_k, m))

    res = scipy_minimize(fun, np.zeros(ncoeff + 1), args=(p_y, cos_k, sin_k, m + 1), method='SLSQP')#, jac=jac)
    print(res.message)
    print("fitting sum of squares error: ", fun(res.x, p_y, cos_k, sin_k, m + 1))#/refs)
    radius = res.x[0]
    coefficients = res.x[1:]
    coefficients = coefficients.reshape((m,2), order = 'F')

    hole.midpoint = midpoint
    return Stellar(midpoint, radius, coefficients)


def fit_stellar(hole, midpoint, refs, angle = 0, ncoeff = 20):
    """

    :param points:
    :param midpoint: new midpoint
    :param N: odd number, number of coeffients for stellar polygon
    :return:
    """
    assert ncoeff % 2 == 0
    # assure points forms a possible stellar polygon with midpoint in midpoint
    # TODO: does not work: if new midpoint is not inside doesnt throw error
    points = hole.discretize_hole(refs)
    points = np.array([np.array(p) for p in points])
    rotated_points = points - np.array(midpoint)
    #if angle != 0:
   # rot_m = rot_matrix(angle)
    #rotated_points = (rot_m @ rotated_points.transpose()).transpose()
    n_points = points.shape[0]
    angles = [math.atan2(p[1],p[0]) for p in rotated_points]
    

    m = int(ncoeff/ 2)

    cos_k = np.array([[np.cos(j * phi)/(j**0.8 + 1) for j in range(1, m+1)] for phi in angles])
    sin_k = np.array([[np.sin(j * phi)/(j**0.8 + 1) for j in range(1, m+1)] for phi in angles])

    p_y = np.array([np.sqrt(rotated_points[i][0] ** 2 + rotated_points[i][1] ** 2) for i in range(n_points)])

   # print(jac(np.zeros(N + 1), p_y, cos_k, sin_k, m))

    res = scipy_minimize(fun, np.zeros(ncoeff + 1), args=(p_y, cos_k, sin_k, m + 1), method='SLSQP')#, jac=jac)
    print(res.message)
    print("fitting sum of squares error: ", fun(res.x, p_y, cos_k, sin_k, m + 1))#/refs)
    radius = res.x[0]
    coefficients = res.x[1:]
    coefficients = coefficients.reshape((m,2), order = 'F')
    print("radius", radius)
    return Stellar(midpoint, radius, coefficients)


def fun(x, *args):
    p_y, cos_k, sin_k, m = args
    n_points = p_y.shape[0]
    f_y = x[0] * np.exp(cos_k @ x[1:m] + sin_k @ x[m:])
    return np.sum((f_y - p_y) ** 2)


def jac(x, *args):
    p_y, cos_k, sin_k, m = args
    n_points = p_y.shape[0]
    a = np.exp(np.ones(n_points) * x[1] + cos_k @ x[2:m + 1] + sin_k @ x[m + 1:])
    d_x_0 = 2 * np.dot(x[0] * a - p_y, a)
    d_x_1 = 2 * x[0] * (np.multiply(cos_k, a[:, np.newaxis])).transpose() @ (x[0] * a - p_y)
    d_x_2 = 2 * x[0] * (np.multiply(sin_k, a[:, np.newaxis])).transpose() @ (x[0] * a - p_y)
    d_x = np.concatenate((np.array([d_x_0]), d_x_1, d_x_2))
    #print(d_x)
    return d_x

def check_valid_midpoint(stellar, npoints = 30):
    midpoint = stellar.midpoint
    points = stellar.discretize_hole(npoints)
    boundary_segments = [LineString([points[i], points[i+1]]) for i in range(npoints-1)]
    boundary_segments.append(LineString([points[-1], points[0]]))
    star_segments = [LineString([points[i], midpoint]) for i in range(npoints)]
    for s in star_segments:
        for b in boundary_segments:
            if b.intersects(s):
                point = b.intersection(s)
                if point.coords[0] in list(s.coords):
                    pass
                else:
                    return False
    return True
