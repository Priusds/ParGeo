from abc import ABC, abstractmethod
import numpy as np
from scipy import optimize
from scipy.optimize import minimize as scipy_minimize
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from ..two_d.star_geometry import *
from ..geo_utils.writeGeo import WriteGeo
from ..geo_utils.utils import write_geo
# TODO: optimize discretization using curve integrals
perm = np.zeros((3,3))
perm[0,0] = 1
perm[1,2] = 1
perm[2,1] = 1

class Bubble(ABC):


    @abstractmethod
    def trafo(self, p):
        """
          Implements the underlying geometric mapping from
          
           p = (theta, phi) in (0,pi) x (0, 2pi) --> F(p) in IR^3

        """
        pass

    @abstractmethod
    def partial_theta(self, p, axis = 'z'):
        pass

    @abstractmethod
    def getPointsOnHeight(self, N, h):
        pass

    #@abstractmethod
    def discretize(self): 
        pass


class Ellipsoid3D(Bubble):


    def __init__(self, midpoint, radii, angles = np.zeros(3)):

        self.M = midpoint # shape 3,1321
        self.a, self.b, self.c = radii
        self.discretize_version = None


        self.alpha, self.beta, self.gamma = angles

        self.Rx  = np.array([ [1,                 0 ,                   0], \
                              [0, np.cos(self.alpha), -np.sin(self.alpha)], \
                              [0, np.sin(self.alpha),  np.cos(self.alpha)]  \
                            ])
        self.Ry  = np.array([ [ np.cos(self.beta), 0, np.sin(self.beta)], \
                              [                 0, 1,                 0], \
                              [-np.sin(self.beta), 0, np.cos(self.beta)]  \
                            ])
        self.Rz  = np.array([ [ np.cos(self.gamma), -np.sin(self.gamma), 0], \
                              [ np.sin(self.gamma),  np.cos(self.gamma), 0], \
                              [                  0,                   0, 1]  \
                            ])
       

    
    def trafo(self, p):
        """
            @param p = (theta, phi)   theta in (0,pi)
                                      phi   in (0,2pi)


                  p[0] is an array of thetas etc

                  p in IR^(2 x N)

        """
        res = np.zeros(3)

        res[0] = (1+0.2*np.sin(3*p[0]) + .3*np.cos(5*p[1]))*self.a * np.sin(p[0]) * np.cos(p[1]) # (1+0.2*np.sin(3*p[0]) + .3*np.cos(5*p[1]))*
        res[1] = self.b * np.sin(p[0]) * np.sin(p[1])
        res[2] = self.c * np.cos(p[0])

        #res  =  self.Rz.dot(self.Ry.dot(self.Rx.dot(res)))

        #res = perm @ res
        res += self.M


        return res

    def partial_theta(self, p, axis = 2):
        """
            Computes the partial derivates w.r.t. theta of the underlying mapping in axis direction

            F(p) = trafo(p) \in  IR^3

            axis == 0 :    delta_theta [F_1](p)
            axis == 1 :    delta_theta [F_2](p)
            axis == 2 :    delta_theta [F_3](p)

        """
        raise NotImplementedError("no partial derivative implemented")

    def getPointsOnHeight2(self, N, h, axis = 2):
        def f(p):
            return self.trafo(p)[axis]

        initial_guesses_theta = np.random.uniform(0,np.pi,N)
        initial_guesses_phi = np.random.uniform(0,2*np.pi,N)
        results = []
       
        for i in range(N):
            sol = optimize.root(f,x0 = [initial_guesses_theta[i], initial_guesses_phi[i]])
            results.append(sol.x)

    
    def getPointsOnHeight(self, N, h, axis = 2):
        """
            Computes points p_i = (theta_i, phi_i), i = 1, ..., N
            
            s.t. trafo(p_i)[axis] = h 


        """

        # check if the height h in IR  intersets the topology:

        #
        #
        #   TODO 
        #
        #


        # create random phi in (0,2pi)    
        #phis = np.random.random(N) * 2*np.pi
        #thetas = np.random.rand(N) * np.pi


        phis = np.linspace(0, 2*np.pi, N)
        thetas = np.linspace(0, np.pi, N)

        def ftheta(phi):
            def f(theta):
                p = np.array([theta, phi])
                return self.trafo(p)[axis] - h
            return f
        def fphi(theta):
            def f(phi):
                p = np.array([theta, phi])
                return self.trafo(p)[axis] - h
            return f

        

        validThetas = []
        validPhis = []

        for phi in phis : 

            f = ftheta(phi)

            Is = np.linspace(0,np.pi,50)
            V = [ f(t) for t in Is]
            minV = min(V)
            maxV = max(V)
            found = False
            I = []
            if min(V) < 0 and max (V) > 0 :
                found = True
                I = Is[V.index(minV)], Is[V.index(maxV)]
                I = [min(I), max(I)]
            if not found:
                #print("Not found")
                continue 
            else:
                #print("FOUND : ", I)
                sol = optimize.root_scalar(f, bracket=I, method='brentq')
                #print(sol.root, sol.iterations, sol.function_calls)
                validPhis.append(phi)
                validThetas.append(sol.root)
        
        for theta in thetas : 

            f = fphi(theta)

            Is = np.linspace(0,2*np.pi,50)
            V = [ f(t) for t in Is]
            minV = min(V)
            maxV = max(V)
            found = False
            I = []
            if min(V) < 0 and max (V) > 0 :
                found = True
                I = Is[V.index(minV)], Is[V.index(maxV)]
                I = [min(I), max(I)]
            if not found:
                #print("Not found")
                continue 
            else:
                #print("FOUND : ", I)
                sol = optimize.root_scalar(f, bracket=I, method='brentq')
                #print(sol.root, sol.iterations, sol.function_calls)
                validThetas.append(theta)
                validPhis.append(sol.root)





        return validThetas, validPhis

    def getHeightRange(self, axis = 2):

        def f(p):
            return self.trafo(p)[axis]
        
        
        resMin = optimize.minimize(f, x0 = np.zeros(2))

        resMax = optimize.minimize(lambda p : -f(p), x0 = np.array([np.pi, 0]))

        hmax, hmin = resMin.fun, (-1) * resMax.fun
        
        if hmin < hmax : 
            return hmin, hmax, resMin.x, resMax.x
        else : 
            return hmax, hmin, resMax.x, resMin.x
        #return min(hmin,hmax), max(hmin,hmax)


    def _cross_section(self, height, nfit):
        """
        finds nfit points at given height and then fits a steller.
        Return fitted stellar.
        """
        points = self.getPointsOnHeight(nfit, height)
        if len(points[0]) < nfit:
            print("found only {0} points instead of {1}".format(len(points), nfit))
            if len(points[0]) < 5:
                raise ValueError("Did not find enough points on height ", height)
        npoints = len(points[0])
        trafo_points = [self.trafo([points[0][i], points[1][i]]) for i in range(npoints)]
        layer_xy = np.array([(p[0], p[1]) for p in trafo_points])
        
        midpoint = get_mass_midpoint(layer_xy)
        print("midpoint", midpoint)
        stellar = fit_stellar(layer_xy, midpoint,ncoeff = 16)
        return stellar

    def discretize_cut_up(self, cut_height, point_id=1, line_id=1,lineLoop_id=1):
        h_min, h_max, p_max, p_min = self.getHeightRange(1)
        print(h_min, h_max)
        assert h_min < cut_height <h_max
        nlayers = 15
        heights = np.linspace(h_min, cut_height,nlayers)
        n = 30
        mid_stellars = []

        for h in heights[1:-2]:
            stellar = self._cross_section(h, n)
            mid_stellars.append(stellar)

        discretize_seq = [40]*(nlayers-2) #[3,4,5,6,5,4,4,3]
        Layers = [s.discretize_hole(d) for s,d in zip(mid_stellars, discretize_seq)]

        Layers = [[(p[0], h, p[1]) for p in l] for l,h in zip(Layers,heights[1:-2])]

        Layers = [[self.trafo(p_min)]] + Layers 

        points, lines, triangles = connect_layers(Layers,point_id,line_id,lineLoop_id)

        self.discretize_version = {"points": points, "lines": lines, "triangles":triangles}
        return self.discretize_version
        
    def discretize_full(self,point_id=1,line_id=1,lineLoop_id=1):
        print("KOOOOOOOOOOOoooo")
        h_min, h_max, p_max, p_min = self.getHeightRange(2)
        nlayers = 24
        heights = np.linspace(h_min, h_max,nlayers)
        n = 40
        mid_stellars = []

        for h in heights[1:-1]:
            stellar = self._cross_section(h, n)
            mid_stellars.append(stellar)

        discretize_seq = [40]*(nlayers-2) #[3,4,5,6,5,4,4,3]
        Layers = [s.discretize_hole(d) for s,d in zip(mid_stellars, discretize_seq)]

        Layers = [[(p[0], h, p[1]) for p in l] for l,h in zip(Layers,heights[1:-1])]
        print("---------------")
        print(len(Layers))
        Layers = [[perm @ self.trafo(p_min)]] + Layers + [[perm @ self.trafo(p_max)]]

        points, lines, triangles = connect_layers(Layers,point_id,line_id,lineLoop_id)

        self.discretize_version = {"points": points, "lines": lines, "triangles":triangles}
        return self.discretize_version
    
    def show(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        self.discretize_full()
        if self.discretize_version is not None:
            for point in self.discretize_version["points"]:
                xs, ys, zs = [], [], []
                for X in point:

                    xs.append(point["x"])
                    ys.append(point["y"])
                    zs.append(point["z"])

                ax.scatter(xs, ys, zs) #, c = 'k') #, ms = 1)


            ax.set_aspect('equal')
            plt.show()
    
    def to_geo_old(self):
        data = self.discretize_full()
        points, lines, triangles = data["points"], data["lines"], data["triangles"]
        file = WriteGeo("test_stellar3d")
        lc = .5
        for point in points:
            file.write_point3d(point["x"], point["y"], point["z"], lc)

        for line in lines:
            file.write_line(line["start"],line["end"])

        for triangle in triangles:
            l1 = triangle["l1"]
            l2 = triangle["l2"]
            l3 = triangle["l3"]
            file.write_line_loop([l1, l2, l3])
            file.write_plane_surface([file.lineLoopId-1])
        file.close()


    def get_geo(self,point_id=1,line_id=1,lineLoop_id=1):
        data = self.discretize_full(point_id,line_id,lineLoop_id)
        points = data["points"]

        lines = data["lines"]
        triangles = data["triangles"]
        p = []
        ll = []
        for i in points:
            Id, x,y,z = i["id"],i["x"],i["y"],i["z"]
            p.append({"id":Id, "coords":[x,y,z],"lc":None})
        for i in triangles:
            Id, l1, l2, l3 = i["id"],i["l1"],i["l2"],i["l3"]
            ll.append({"id":Id, "lines":[l1,l2,l3]})

        surfaces = [{"id":l["id"],"lineLoops":[l["id"]]} for l in ll]
        geo = {"points":p,"lines":lines,"lineLoops":ll,"surfaces":surfaces,"surfaceLoops":[],"volumes":[]}
        return geo

    def to_geo(self, file_name,point_id=1,line_id=1,lineLoop_id=1):
        data = self.discretize_full(point_id,line_id,lineLoop_id)
        points = data["points"]

        lines = data["lines"]
        triangles = data["triangles"]
        p = []
        ll = []
        for i in points:
            Id, x,y,z = i["id"],i["x"],i["y"],i["z"]
            p.append({"id":Id, "coords":[x,y,z],"lc":None})
        for i in triangles:
            Id, l1, l2, l3 = i["id"],i["l1"],i["l2"],i["l3"]
            ll.append({"id":Id, "lines":[l1,l2,l3]})

        surfaces = [{"id":l["id"],"lineLoops":[l["id"]]} for l in ll]
        geo = {"points":p,"lines":lines,"lineLoops":ll,"surfaces":surfaces,"surfaceLoops":[],"volumes":[]}
        write_geo(geo, file_name)

def locate(geo, lpoint_id,h=-.06):
    # TODO: point ID not -1
    #self.geo["point_id"], self.geo["line_id"], self.geo["lineLoop_id"]print(geo["points"])
    lpoint_index = [search_index("id",Id,geo["points"]) for Id in lpoint_id]
    lpoint = [geo["points"][i] for i in lpoint_index]
    loc = [p["coords"][1] < h for p in lpoint] # false if out
    if sum(loc) == 3:
        return 1,loc
    elif sum(loc) == 0:
        return -1,loc
    else:
        return 0,loc

def locate_point(geo, point, h):
    point_index = search_index("id", point,geo["points"])
    return geo["points"][point_index]["coords"][1] < h

def sort_lineLoop(geo, lines_id):
    lines = [geo["lines"][search_index("id",l,geo["lines"])] for l in lines_id]

    #lines = [geo["lines"][l-1] for l in lines_id]
    nlines = len(lines)
    free_lines = lines_id
    sorted_lines = [free_lines[0]]
    free_lines.pop(0)

    i = 0
    geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["start"] #geo["lines"][sorted_lines[-1]-1]["start"]

    while len(sorted_lines) < nlines and i < 1000:
        i+=1
        current = set([geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["start"],
                        geo["lines"][search_index("id", sorted_lines[-1], geo["lines"]) ]["end"]])
        #current = set([geo["lines"][sorted_lines[-1]-1]["start"],geo["lines"][sorted_lines[-1]-1]["end"]])

        cond = True
        j = 0
        
        while cond and j < len(free_lines):
            free_line = set([geo["lines"][search_index("id", free_lines[j], geo["lines"])]["start"], geo["lines"][search_index("id", free_lines[j], geo["lines"])]["end"]])
            #free_line = set([geo["lines"][free_lines[j]-1]["start"], geo["lines"][free_lines[j]-1]["end"]])
            if len(free_line.intersection(current)) > 0:

                sorted_lines.append(free_lines[j])
                free_lines.remove(free_lines[j])
                cond = False
            j += 1

    if len(sorted_lines) != nlines:
        raise ValueError("Problem with clipping: no sorting for lineLoop")


    for i in range(1,nlines):
        current_start = geo["lines"][search_index("id",abs(sorted_lines[i]), geo["lines"])]["start"]
        #current_start = geo["lines"][sorted_lines[i]-1]["start"]
        if sorted_lines[i-1] < 0:
            prev_end = geo["lines"][search_index("id",abs(sorted_lines[i-1]), geo["lines"])]["start"]
            #prev_end = geo["lines"][abs(sorted_lines[i-1])-1]["start"]
        else:
            prev_end = geo["lines"][search_index("id",abs(sorted_lines[i-1]), geo["lines"])]["end"]
            #prev_end = geo["lines"][abs(sorted_lines[i-1])-1]["end"]
        if prev_end != current_start:
            sorted_lines[i] = -sorted_lines[i]
        

    return sorted_lines

def search_index(key, value, list_of_dicts):
    l = len(list_of_dicts)
    i = 0
    while i < l:
        if list_of_dicts[i][key] == value:
            return i
        i += 1
    raise ValueError("Did not find ", value, " of ", key)

def clip(geo, h):
    finalPoints = []
    finalLines = []
    finalLineLoops = []
    clip_geo = {"points":finalPoints, "lines":finalLines, "lineLoops":finalLineLoops}
    triangles_total = geo["lineLoops"]
    triangles_inside = []
    triangles_outside = []
    triangles_inout = []
    L_out = []
    L_on = []
    L_in = []
    P_out = []
    P_on = []
    P_in = []
    P_ = []
    lineLoop = []
    i = 1

    for triangle in geo["lineLoops"]:
        
        l1, l2, l3 = triangle["lines"]
        lines = [geo["lines"][search_index("id",abs(l),geo["lines"])] for l in [l1,l2,l3]]
        #print([geo["lines"][l-1] for l in [l1,l2,l3]])
        #points = list(set([geo["lines"][abs(l)-1]["start"] for l in triangle["lines"]] + [geo["lines"][abs(l)-1]["end"] for l in triangle["lines"]]))
        points = list(set([l["start"] for l in lines] + [l["end"] for l in lines]))
        
        loc_sum,loc = locate(geo, points, h) # 1 in, 0 mix, -1 out
        
        if loc_sum == 1: # IN
            triangles_inside.append(triangle["id"])
            L_in += [abs(l) for l in triangle["lines"]]
            P_in += points
        elif loc_sum == -1: # OUT
            triangles_outside.append(triangle["id"])
            P_out += points
            L_out += [abs(l) for l in triangle["lines"]]
        else: # IN/OUT
            lines_on = [abs(l) for l in triangle["lines"]]
            L_on += lines_on
            triangles_inout.append(triangle["id"])
            P_on += points
            for i,boolean in enumerate(loc):
                if not boolean:
                    P_.append(points[i])
            
            for l in lines_on:
                p1 = geo["lines"][search_index("id",l,geo["lines"])]["start"]
                p2 = geo["lines"][search_index("id",l,geo["lines"])]["end"]
                #p1 = geo["lines"][l-1]["start"]
                #p2 = geo["lines"][l-1]["end"]
                if ((not locate_point(geo, p1, h)) and (not locate_point(geo, p2, h))):
                    lineLoop.append(l)


    for triangle in triangles_inout:
        triangle_index = search_index("id",triangle,geo["lineLoops"])
        lines = geo["lineLoops"][triangle_index]["lines"] # lines = geo["lineLoops"][triangle -1 ]["lines"] 
        lines = [abs(l) for l in lines]
        L_in += lines

    L_out = set(L_out)
    L_on = set(L_on)
    L_in = set(L_in)
    P_out = set(P_out)
    P_in = set(P_in)
    P_on = set(P_on)

    L_tot = L_in # .union(L_out.intersection(L_on))


    finalPoints += [geo["points"][search_index("id",i,geo["points"])] for i in list(P_in.union(P_on))]
    #finalPoints += [geo["points"][i-1] for i in list(P_in.union(P_on))]

    # Process finalPoints
    for p in P_:
        for pp in finalPoints:
            if pp["id"] == p:
                coords = pp["coords"]
                pp["coords"] = [coords[0], h ,coords[2]]
    
    finalLines += [geo["lines"][search_index("id",i,geo["lines"])] for i in list(L_tot)]
    #finalLines += [geo["lines"][i-1] for i in list(L_tot)]
    finalLineLoops += [geo["lineLoops"][search_index("id",i,geo["lineLoops"])] for i in (triangles_inside + triangles_inout)]
    #finalLineLoops += [geo["lineLoops"][i-1] for i in (triangles_inside + triangles_inout)]

    # Sort line Loop

    
    sorted_lineLoop = sort_lineLoop(geo, list(set(lineLoop)))
    
    intersectionLoop = {"id":len(geo["lineLoops"])+1,"lines":sorted_lineLoop}
    #clip_geo["lineLoops"].append(intersectionLoop)
    #clip_geo["surfaces"] = [{"id":ll["id"], "lineLoops":ll["id"]} for ll in clip_geo["lineLoops"]]
    #write_geo(clip_geo, "test_bubbly")
    return clip_geo, intersectionLoop  


def visualizeLayers(stellar3D, N, hrange, noLayers, axis = 2):

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        #H = [-0.99,-0.9,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.750, 0.9,0.99]

        H = np.linspace(hrange[0], hrange[1], noLayers)

        #H = H[1:-2]

        for h in H:
            thetas, phis = stellar3D.getPointsOnHeight(N, h, axis)

            xs, ys, zs = [], [], []

            for i in range(len(phis)):
                p = np.array([thetas[i], phis[i]])
                #p = np.concatenate(([thetas[i]],[phis[i]]))
                X = stellar3D.trafo(p)

                xs.append(X[0])
                ys.append(X[1])
                zs.append(X[2])
            print(zs)
            ax.scatter(xs, ys, zs, c = 'k')#, ms = 1)


        #ax.set_aspect('equal')
        plt.show()

def fun(x, *args):
    p_y, cos_k, sin_k, m = args
    n_points = p_y.shape[0]
    f_y = x[0] * np.exp(cos_k @ x[1:m] + sin_k @ x[m:])
    return np.sum((f_y - p_y) ** 2)

def fit_stellar(points, midpoint, angle = 0, ncoeff = 20):
    """

    :param points: numpy array of shape (2,N) TODO: check
    :param midpoint: new midpoint
    :param N: odd number, number of coeffients for stellar polygon
    :return:
    """
    assert ncoeff % 2 == 0
    # assure points forms a possible stellar polygon with midpoint in midpoint
    # TODO: does not work: if new midpoint is not inside doesnt throw error
    #print(midpoint[:,0].shape)
    rotated_points = points - midpoint
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

def get_mass_midpoint(points):
    assert points.shape[1] == 2
    return np.mean(points, 0)

def get_stellar(points):
    midpoint = get_mass_midpoint(points)
    stellar = fit_stellar(points, midpoint)
    if not stellar.check_valid():
        raise ValueError("Not possible to fit stellar")
    else:
        return stellar

def connect_layers(layers,current_point_ID, current_line_ID,current_triangle_ID):
    # layer is a list of lists
    # every sublist contains the points on a certain hight
    points = [] # point = {"id", "x", "y", "z"}
    lines = []  # line = {"id", "start", "end"}
    triangles = [] # triangle {"id", "l1","l2","l3"}
    nlayers = len(layers)
    lnpoints = [len(l) for l in layers]
    start_point_ID = current_point_ID
    start_line_ID = current_line_ID
    start_triangle_ID = current_triangle_ID

    temp = np.cumsum(lnpoints)
    lpoint_ID = temp + current_point_ID
    lpoint_ID = [current_point_ID] + list(lpoint_ID)
    lline_ID = temp + current_line_ID - 1
    lline_ID = list(lline_ID)

    # Add all points
    for i in range(nlayers):
        for j in range(lnpoints[i]):
            x,y,z = layers[i][j]
            point = {"id":current_point_ID,"x":x, "y":y,"z":z}
            points.append(point)
            current_point_ID += 1
    
    for i in range(nlayers):
        for j in range(lnpoints[i]):
            if lnpoints[i] >= 3:
                t = sum([lnpoints[k] for k in range(i)])
                s = (j+1)%lnpoints[i]
                line = {"id":current_line_ID, "start":start_point_ID+t+j,"end": start_point_ID+t+s}
                lines.append(line)
                current_line_ID += 1

    for i in range(nlayers-1):
        start_point_ID_1 = lpoint_ID[i]
        start_point_ID_2 = lpoint_ID[i+1]
        start_point_ID_3 = lpoint_ID[i+2]
        

        lines.append({"id":current_line_ID, "start":start_point_ID_1,"end":start_point_ID_2})

        current_line_ID += 1
        start_line_ID = current_line_ID - 1
        next_point_ID_1 = start_point_ID_1 + 1
        next_point_ID_2 = start_point_ID_2 + 1
        choose_triangle = True
        better_triangle = None
        iter = 0
       # print("start: ",start_point_ID_1, start_point_ID_2, start_point_ID_3)
        while ((next_point_ID_1 < start_point_ID_2) or (next_point_ID_2 < start_point_ID_3)) and iter <= 1000:
            iter += 1
            if choose_triangle:
                better_triangle = np.random.choice([1,2])
                if (better_triangle == 1 or next_point_ID_2 >= start_point_ID_3) and next_point_ID_1 < start_point_ID_2:
                    line = {"id": current_line_ID, "start": next_point_ID_1,"end": next_point_ID_2 - 1}
                    lines.append(line)
                    current_line_ID += 1
                    triangle = {"id":current_triangle_ID,"l1":current_line_ID-2,"l2":-(current_line_ID-1),"l3":-(lline_ID[i-1]+next_point_ID_1 - start_point_ID_1-1)}
                    triangles.append(triangle)
                    current_triangle_ID += 1
                    next_point_ID_1 += 1
                else:
                    line = {"id": current_line_ID, "start": next_point_ID_1 - 1,"end": next_point_ID_2 }
                    lines.append(line)
                    current_line_ID += 1
                    triangle = {"id":current_triangle_ID,"l1":current_line_ID-1,"l2":-(lline_ID[i] + next_point_ID_2 - start_point_ID_2-1),"l3":-(current_line_ID-2)}
                    triangles.append(triangle)
                    current_triangle_ID += 1
                    next_point_ID_2 += 1
                    better_triangle = 2
                choose_triangle = False
            else:

                if better_triangle == 1:
                    if next_point_ID_2 < start_point_ID_3:
                        line = {"id": current_line_ID, "start":next_point_ID_1-1,"end":next_point_ID_2}
                        lines.append(line)
                        current_line_ID += 1
                        triangle = {"id":current_triangle_ID,
                                    "l1":current_line_ID-1,
                                    "l2":-(lline_ID[i] + next_point_ID_2 - start_point_ID_2-1),
                                    "l3":-(current_line_ID-2)}
                        triangles.append(triangle)
                        current_triangle_ID += 1
                        next_point_ID_2 += 1
                else:
                    if next_point_ID_1 < start_point_ID_2:
                        line = {"id": current_line_ID, "start":next_point_ID_1,"end":next_point_ID_2-1}
                        lines.append(line)
                        current_line_ID += 1
                        triangle = {"id":current_triangle_ID,
                                    "l1":current_line_ID-2,
                                    "l2":-(current_line_ID-1),
                                    "l3":-(lline_ID[i-1]+next_point_ID_1 - start_point_ID_1-1)}
                        triangles.append(triangle)
                        current_triangle_ID += 1
                        next_point_ID_1 += 1
                choose_triangle = True

    
        if next_point_ID_1 == start_point_ID_2 and next_point_ID_2 == start_point_ID_3:
            # TODO: improve triangle selection, angle based
            better_triangle =np.random.choice([1,2])
            if lnpoints[i] == 1:
                better_triangle = 2
                triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":(lline_ID[i] + next_point_ID_2 - start_point_ID_2-1),
                            "l3":-start_line_ID}
                triangles.append(triangle)
                current_triangle_ID += 1
            elif lnpoints[i+1] == 1:
                better_triangle = 1
                triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":-start_line_ID,
                            "l3":-(lline_ID[i]-1)} #-(lline_ID[i] + next_point_ID_2 - start_point_ID_2-1)}
                triangles.append(triangle)
                current_triangle_ID += 1

            else:
                if better_triangle == 1:
                    line = {"id":current_line_ID,"start":start_point_ID_1,"end":start_point_ID_3-1}
                    lines.append(line)
                    current_line_ID += 1

                    triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":-(current_line_ID-2),
                            "l3":(lline_ID[i]-1)} 
                    triangles.append(triangle)
                    current_triangle_ID += 1
                    
                    triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":(lline_ID[i+1]-1),
                            "l3":-start_line_ID} 
                    triangles.append(triangle)
                    current_triangle_ID += 1

                    
                else:
                    line = {"id":current_line_ID,"start":start_point_ID_2-1,"end":start_point_ID_2}
                    lines.append(line)
                    current_line_ID += 1

                    triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":-(lline_ID[i+1]-1),
                            "l3":-(current_line_ID-2)} 
                    triangles.append(triangle)
                    current_triangle_ID += 1


                    triangle = {"id":current_triangle_ID,
                            "l1":current_line_ID-1,
                            "l2":-start_line_ID,
                            "l3":-(lline_ID[i]-1)} 
                    triangles.append(triangle)
                    current_triangle_ID += 1

    return points, lines, triangles 

    