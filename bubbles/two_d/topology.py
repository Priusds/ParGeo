"""Core module for the creation of 2d geometries."""
import json
from typing import Optional

from bubbles.geo_utils.utils import write_geo
from bubbles.two_d.hole import Circle, Ellipse, Stellar
from bubbles.two_d.intersections import (
    convex_hull,
    intersection_hole_corner,
    intersection_hole_edge,
)
from bubbles.two_d.sat import collide

# ---------- Version 2 -------------


class Topology(object):
    """Main class for the creation of two-dimensional perforated geometries.

    TODO:
    Steps for the refactoring of this class:
        - Understand present functionalities:
            - What are minimal topology assumptions, s.t. the code
                works Bug-free?
            - What is the coolest thing we can do, but only works sometimes?

        - Add visualization functionalities, based on matplotlib:
            - Given any type of geometry object (actually the class Hole
                could be generalized to a 2d-geometry class, if handy maybe
                even subclass some geometry package) write a plot functionality
                that always plots some current topology (super easy).
                For debugging this is really needed.
        
        - State the foundamental needed output of this class / bubbles.two_d package:
            - Create objects that 
    
    Class description:
    A geometry is defined by following elements:
        - domain
            - Has one line-loop
            - Has one Plane Surface and one Physical Surface
        - sub-domains
            - Has two line-loops (inner and outer rect)
            - Has one Plane Surface and one Physical Surface
        - holes (aka bubbles)
            - A filled holes gets his own line loop, if it intersects with domain edges.
            - How is a domain filled / empty in free space? (Just by changing the orientation?) # TODO: check this
            - If filled: has it's own Plane Surface
            - If intersects with subdomain: split into two "different" holes, inner and outer one.
            
    For now domain and sub-domains are rectangles (different shapes might be 
    allowed in future). The holes are objects from sub-classes of `bubble.two-d.hole.Hole`.

    Workflow:
    1.) Instantiate a topology object => defines the domain and sub-domains.
    2.) Add holes via self.add_hole(...) => adds a hole to the topology.
        2.1) First it's checked wether the hole can be added. Criterias are:
            - Holes do not intersect other holes.
            - Holes cannot intersect opposite domain (also sub-domain?) edges.
            - Holes cannot intersect adjactent edges without not-containing the corner-vertex.
            - (Basically the domain must remain connected. Same logic applies for sub-domains?)
            - TODO: Check what is missing / correct.
    3.) Write .GEO file via self.get_geometry(...) and `bubbles.geo_utils.utils.write_geo`
        3.1) The method self.get_geometry(...) creates all the linking between
        point-ID/line-ID/ect.
    """
    def __init__(self, rect_out, lrect_in, periodic_boundary=True, method="linear"):
        # parameters for rect_out
        x0_out, y0_out = rect_out[0]
        x1_out, y1_out = rect_out[1]
        self.x0_out = x0_out
        self.y0_out = y0_out
        self.x1_out = x1_out
        self.y1_out = y1_out

        self.dx = x1_out - x0_out
        self.dy = y1_out - y0_out

        self.nRects = len(lrect_in)
        # parameters for rect_in
        self.periodic_boundary = periodic_boundary

        # underlying method of hole and rectangle intersection point computation
        # method = 'linear' / 'optimise'
        # linear   : computes intersection points based in segment-edge intersection based on holes discretization points
        # optimise : computes apprx to exact intersection of hole and edge based on numerical approximation
        self.method = method

        self.rects_in = {
            Id
            + 1: {
                "x0": lrect_in[Id][0][0],
                "y0": lrect_in[Id][0][1],
                "x1": lrect_in[Id][1][0],
                "y1": lrect_in[Id][1][1],
                "holes_in": [],
                "edgeLeft": [
                    [(lrect_in[Id][0][0], lrect_in[Id][1][1])]
                ],  # This is a list of holes
                "edgeUp": [
                    [(lrect_in[Id][1][0], lrect_in[Id][1][1])]
                ],  # each hole can have either length 1 or 3
                "edgeRight": [
                    [(lrect_in[Id][1][0], lrect_in[Id][0][1])]
                ],  # length == 1 means: vertex case
                "edgeDown": [
                    [(lrect_in[Id][0][0], lrect_in[Id][0][1])]
                ],  # length == 3 means: [intersection_points, inner_points, outer_points] (ccw)
            }
            for Id in range(len(lrect_in))
        }
        self.rects_in[0] = {
            "x0": x0_out,
            "y0": y0_out,
            "x1": x1_out,
            "y1": y1_out,
            "edgeLeft": [[(x0_out, y1_out)]],
            "edgeUp": [[(x1_out, y1_out)]],
            "edgeRight": [[(x1_out, y0_out)]],
            "edgeDown": [[(x0_out, y0_out)]],
        }

        # holes between rect_in and rect_out (no intersections)
        self.holes_between = []

        # holes on rect_out (intersection with rect_out)
        self.edgeLeft_out = [[(x0_out, y1_out)]]
        self.edgeUp_out = [[(x1_out, y1_out)]]
        self.edgeRight_out = [[(x1_out, y0_out)]]
        self.edgeDown_out = [[(x0_out, y0_out)]]

        # list of all holes (convex hull of the hole) [checking for intersections]

        self.all_holes_conv_hull = []

        # description of the topology in json format (more readable)

        self.topology_json = {
            "rect_out": rect_out,
            "rects_in": lrect_in,
            "boundaryCondition": periodic_boundary,
        }

    def add_hole(self, hole, refs=60):
        # TODO: Missing corner case
        assert (
            isinstance(hole, Circle)
            or isinstance(hole, Ellipse)
            or isinstance(hole, Stellar)
        )
        discretized_hole = hole.discretize_hole(refs)

        if self.__add_hole((discretized_hole, hole)):
            counter = len(self.topology_json) - 3
            self.topology_json["hole_" + str(counter)] = (hole.to_dict(), refs)
            print("Hole added")
            return True
        print("Hole could not be added")
        return False

    def write_geo(self, file_name, filled=True, lc_in=0.2, lc_out=0.2):
        """writes a .geo file of the current topology"""
        geometry = self.get_geometry(lc_in, lc_out, filled)
        write_geo(geometry, file_name)

    def safe_json(self, file_name):
        # raise NotImplementedError("Change from 1 inner rect to multiple is missing")
        """safes the current topology as json file (without refinements)"""
        with open(file_name + ".json", "w") as fp:
            json.dump(self.topology_json, fp, indent=2)

    @staticmethod
    def load_topology(file_name, hole_refinement=None, decay=None):
        # raise NotImplementedError("Change from 1 inner rect to multiple is missing")
        """
        refinement_pattern = default_refinement_pattern (constant verfeinerung oder sowas)


        Extracts all necessary infos from topology_json to build topology information

        and via refinement_pattern (that info needs to be stored here also)
        calling writeFile() builds geo-File.


        hole_refinement should be a method with 3 inputs:
            rect_out
            rect_in
            hole_infos (for example midpoint)
        """

        with open(file_name + ".json", "r") as fp:
            topo_json = json.load(fp)

        topo = Topology(
            topo_json["rect_out"], topo_json["rects_in"], topo_json["boundaryCondition"]
        )

        for data in topo_json.items():
            key, val = data
            if "hole" in key:
                hole, refs = val
            # for hole in topo_json.values():

            if type(hole) is dict:
                if hole["type"] == "circle":
                    hole = Circle(hole["midpoint"], hole["radius"])
                    hole_d = hole.discretize_hole(refs)
                    topo.__add_hole(hole_d)
                    # if hole_refinement is None:
                    #    hole_d = hole.discretize_hole(refs)
                    #    topo.__add_hole(hole_d)
                    # else:
                    #    hole_d = hole.discretize_hole(
                    #        hole_refinement(topo_json["rect_out"], topo_json["rect_in"], hole.midpoint, decay))
                    #    topo.__add_hole(hole_d)
                else:
                    raise NotImplementedError("Other hole types need to be specified.")
                # elif hole['type'] == 'ellipse':
                #    hole = Ellipse(hole['midpoint'], hole['axis'], hole['angle'])
                #    if hole_refinement is None:
                #        hole_d = hole.discretize_hole(refs)
                #        topo.__add_hole(hole_d)
                #    else:
                #        hole_d = hole.discretize_hole(
                #            hole_refinement(topo_json["rect_out"], topo_json["rect_in"], hole.midpoint, decay))
                #        topo.__add_hole(hole_d)

                # elif hole['type'] == 'stellar':
                #    hole = Stellar(hole['midpoint'], hole['radius'], hole['coefficient'])
                #    if hole_refinement is None:
                #        hole_d = hole.discretize_hole(4)#+++++++++++++++++++++++++++++++
                #        topo.__add_hole(hole_d)
                #    else:
                #        hole_d = hole.discretize_hole(
                #            hole_refinement(topo_json["rect_out"], topo_json["rect_in"], hole.midpoint, decay))
                #        topo.__add_hole(hole_d)

                # elif hole['type'] == 'rectangular':
                #    hole = Rectangular(hole['midpoint'], hole['width'], hole['height'])
                #    topo.__add_hole(hole.discretize_hole(None))

        return topo

    def __add_hole(self, hole):
        """
        This is the main method of Topology.

        :param hole: list. Each element of hole shoud be a lenght 2 list / tuple
        :return: True if the hole could be added, False if not

        The input hole should consist of the vertices (of the discretized hole) in ccw order.

        example: hole = [(0,0), (1,0), (1,1), (0,1)]

        """
        dhole = hole[0]
        # Check intersection with other holes
        if len(dhole) <= 2:
            # The hole must consist of at least 3 points
            return False
        convex_hull_hole = convex_hull(dhole)  # dhole #
        if not self.__valid_input(convex_hull_hole):
            return False

        # localize the hole and check if the intersections with rect_in and rect_out is regular
        # (more checks will follow)
        # intersection_loc in {None, "left","left_down",...}
        # rect_id in {0,1,2,3...} 0 is rect_out, i > 0 is the ID of one rect_in
        bool__, intersection_loc, rect_id = self.__localize(hole)

        if bool__:
            if intersection_loc is None:
                # ---------------NO INTERSECTIONS-------------------
                if rect_id > 0:
                    # hole is in rect_in
                    self.rects_in[rect_id]["holes_in"].append(dhole)

                else:
                    # hole is between rect_in and rect_out
                    self.holes_between.append(dhole)

            elif intersection_loc in {"left", "up", "right", "down"}:
                # ----------------EDGE INTERSECTION ----------------
                self.__add_hole_on_edge(
                    intersection_loc, rect_id, hole, convex_hull_hole
                )

            else:
                # ---------------CORNER INTERSECTION ------------------
                self.__add_hole_on_corner(intersection_loc, rect_id, hole)

            self.all_holes_conv_hull.append(convex_hull_hole)
            return True
        else:
            print(intersection_loc)
            return False

    def __valid_input(self, new_hole):
        """
        Checks if the convex hull of new_hole intersects with
        the convex hull of another hole
        :param new_hole: list of vertices in ccw order
        :return: True if there aren't any intersections, False in the other case
        """
        for hole in self.all_holes_conv_hull:
            if collide(new_hole, hole):
                return False
        return True

    def __valid_reflection_corner(self, hole_d_l, hole_u_l, hole_u_r, hole_d_r):
        # For now hole_d_l doesn't contain the intersection points, therefore convex_hull of less then 3 points
        # would make an error
        if min([len(hole_d_l), len(hole_u_r), len(hole_u_l), len(hole_d_r)]) <= 2:
            return False

        d_l = convex_hull(hole_d_l)
        u_l = convex_hull(hole_u_l)
        u_r = convex_hull(hole_u_r)
        d_r = convex_hull(hole_d_r)
        d_l.append([self.x1_out, self.y1_out])
        u_l.append([self.x1_out, self.y0_out])
        u_r.append([self.x0_out, self.y0_out])
        d_r.append([self.x0_out, self.y1_out])

        if (
            self.__valid_input(d_l)
            and self.__valid_input(u_l)
            and self.__valid_input(u_r)
            and self.__valid_input(d_r)
        ):
            self.all_holes_conv_hull.append(d_l)
            self.all_holes_conv_hull.append(u_l)
            self.all_holes_conv_hull.append(u_r)
            self.all_holes_conv_hull.append(d_r)
            return True
        return False

    def __localize(self, hole) -> tuple[bool, Optional[str], Optional[int]]:
        """ " localizes the hole
            return bool, str, int
                -bool: if the hole can be added or not
                -str: location of intersection edge/edges str in {"left","left_down",None,...}
                -int: None if no Intersection, 0 if intersection with rect_out, Id if intersection with some rect in
        input: list of tuples
        output: bool_ean, str, str

        Possible outcomes:

            False, "outside domain", None
            True, "right", rect_in_id # or 0
            True, "left", rect_in_id # or 0
            True, "down", rect_in_id # or 0
            True, "up", rect_in_id # or 0
            True, "down_left", rect_in_id # or 0
            True, "left_up", rect_in_id # or 0
            True, "right_down", rect_in_id # or 0
            True, "up_right", rect_in_id # or 0
            False, "Doesnt intersects with a corner", None
            False, "Intersection with too many edges, hole cannot be added", None
            False, "Intersection with rect_in and rect_out", None
            False, "Intersection with multiple rects inside", None
            False, "no proper corner intersection", None
            True, None, rect_in_id
        """

        dhole = hole[0]

        # TODO make it possible to choose how to compute intersection points
        max_x = max(dhole, key=lambda x: x[0])[0]
        max_y = max(dhole, key=lambda x: x[1])[1]

        min_x = min(dhole, key=lambda x: x[0])[0]
        min_y = min(dhole, key=lambda x: x[1])[1]

        x0_out = self.x0_out
        x1_out = self.x1_out
        y0_out = self.y0_out
        y1_out = self.y1_out

        if max_x < x0_out or min_x > x1_out or max_y < y0_out or min_y > y1_out:
            return False, "outside domain", None

        elif min_x < x0_out or max_x > x1_out or min_y < y0_out or max_y > y1_out:
            # ----------------INTERSECTION RECT OUT---------------------
            intersection_with_rects_in = False
            for rect_in_id, rect_in in self.rects_in.items():
                if rect_in_id > 0:
                    x0_in = rect_in["x0"]
                    y0_in = rect_in["y0"]
                    x1_in = rect_in["x1"]
                    y1_in = rect_in["y1"]
                    if max_x < x0_in or min_x > x1_in or max_y < y0_in or min_y > y1_in:
                        # no intersection bewteen hole and rect_in
                        pass
                    else:
                        intersection_with_rects_in = True
            if not intersection_with_rects_in:
                v = (min_x < x0_out, max_x > x1_out, min_y < y0_out, max_y > y1_out)
                # hole doesnt intersects with rect in
                if sum(v) == 1:
                    # ----------- INTERSECTION ONE EDGE ----------------------

                    if v[0] == 1:
                        return True, "left", 0
                    elif v[1] == 1:
                        return True, "right", 0
                    elif v[2] == 1:
                        return True, "down", 0
                    else:
                        return True, "up", 0

                elif sum(v) == 2:
                    # ------------INTERSECTION TWO EDGES---------------------

                    if min_x < x0_out and min_y < y0_out:
                        return True, "down_left", 0

                    elif min_x < x0_out and max_y > y1_out:
                        return True, "left_up", 0

                    elif max_x > x1_out and min_y < y0_out:
                        return True, "right_down", 0

                    elif max_x > x1_out and max_y > y1_out:
                        return True, "up_right", 0

                    else:
                        return False, "Doesnt intersects with a corner", None
                else:
                    return (
                        False,
                        "Intersection with too many edges, hole cannot be added",
                        None,
                    )
            else:
                return False, "Intersection with rect_in and rect_out", None

        else:
            # ---------------- NO INTERSECTION RECT OUT ------------------
            # Now check if hole intersects with inner rects
            # Only allowed are no intersections with rects in or exactly one
            lcollisions = []
            for rect_in_id, rect_in in self.rects_in.items():
                if rect_in_id > 0:
                    x0_in = rect_in["x0"]
                    y0_in = rect_in["y0"]
                    x1_in = rect_in["x1"]
                    y1_in = rect_in["y1"]

                    if max_x < x0_in or min_x > x1_in or max_y < y0_in or min_y > y1_in:
                        # hole does not intersect with rect_in
                        pass
                    elif (
                        min_x < x0_in or max_x > x1_in or min_y < y0_in or max_y > y1_in
                    ):
                        # hole intersect with rect in
                        lcollisions.append((rect_in_id, "on"))
                    else:
                        # hole is inside rect in
                        lcollisions.append((rect_in_id, "in"))

            if len(lcollisions) > 1:
                return False, "Intersection with multiple rects inside", None

            elif len(lcollisions) == 0:
                # hole is inside rect out and does not intersect with inner rects
                return True, None, 0
            else:
                rect_in_id = lcollisions[0][0]
                rect_in = self.rects_in[rect_in_id]
                x0_in = rect_in["x0"]
                y0_in = rect_in["y0"]
                x1_in = rect_in["x1"]
                y1_in = rect_in["y1"]
                loc = lcollisions[0][1]

                if loc == "on":
                    # ------------------INTERSECTION RECT IN------------------

                    v = (min_x < x0_in, max_x > x1_in, min_y < y0_in, max_y > y1_in)
                    if sum(v) == 1:
                        # --------------- INTERSECTION ONE EDGE-----------------

                        if v[0] == 1:
                            return True, "left", rect_in_id
                        elif v[1] == 1:
                            return True, "right", rect_in_id
                        elif v[2] == 1:
                            return True, "down", rect_in_id
                        else:
                            return True, "up", rect_in_id

                    elif sum(v) == 2:
                        # ---------------- INTERSECTION TWO EDGES---------------------

                        if min_x < x0_in and min_y < y0_in:
                            # -----------------DOWN LEFT------------------
                            bool_s, intersection_s = intersection_hole_edge(
                                hole, (x0_in, "x"), self.method
                            )
                            bool_r, intersection_r = intersection_hole_edge(
                                hole, (y0_in, "y"), self.method
                            )
                            if bool_r and bool_s:
                                s0, s1 = intersection_s
                                r0, r1 = intersection_r
                                if s0[1] < y0_in < s1[1] or s0[1] > y0_in > s1[1]:
                                    if r0[0] < x0_in < r1[0] or r0[0] > x0_in > r1[0]:
                                        return True, "down_left", rect_in_id
                                    else:
                                        return (
                                            False,
                                            "no proper corner intersection",
                                            None,
                                        )

                                elif (
                                    max(s0[1], s1[1]) < y0_in
                                    and min(r0[0], r1[0]) > x0_in
                                ):
                                    return True, "down", rect_in_id

                                elif (
                                    max(r0[0], r1[0]) < x0_in
                                    and min(s0[1], s1[1]) > y0_in
                                ):
                                    return True, "left", rect_in_id
                                else:
                                    return False, "no proper corner intersection", None

                            else:
                                return False, "no proper corner intersection", None

                        elif min_x < x0_in and max_y > y1_in:
                            # ---------------LEFT UP---------------------
                            bool_s, intersection_s = intersection_hole_edge(
                                hole, (x0_in, "x"), self.method
                            )
                            bool_r, intersection_r = intersection_hole_edge(
                                hole, (y1_in, "y"), self.method
                            )

                            if bool_r and bool_s:
                                s0, s1 = intersection_s
                                r0, r1 = intersection_r

                                if s0[1] < y1_in < s1[1] or s0[1] > y1_in > s1[1]:
                                    if r0[0] < x0_in < r1[0] or r0[0] > x0_in > r1[0]:
                                        return True, "left_up", rect_in_id
                                    else:
                                        return (
                                            False,
                                            "no proper corner intersection",
                                            None,
                                        )

                                elif (
                                    max(r0[0], r1[0]) < x0_in
                                    and max(s0[1], s1[1]) < y1_in
                                ):
                                    return True, "left", rect_in_id

                                elif (
                                    min(r0[0], r1[0]) > x0_in
                                    and min(s0[1], s1[1]) > y1_in
                                ):
                                    return True, "up", rect_in_id

                                else:
                                    return False, "no proper corner intersection", None

                            else:
                                return False, "no proper corner intersection", None

                        elif max_x > x1_in and min_y < y0_in:
                            # -------------RIGHT DOWN-----------------

                            bool_s, intersection_s = intersection_hole_edge(
                                hole, (x1_in, "x"), self.method
                            )
                            bool_r, intersection_r = intersection_hole_edge(
                                hole, (y0_in, "y"), self.method
                            )
                            if bool_r and bool_s:
                                s0, s1 = intersection_s
                                r0, r1 = intersection_r

                                if s0[1] < y0_in < s1[1] or s0[1] > y0_in > s1[1]:
                                    if r0[0] < x1_in < r1[0] or r0[0] > x1_in > r1[0]:
                                        return True, "right_down", rect_in_id
                                    else:
                                        return (
                                            False,
                                            "no proper corner intersection",
                                            None,
                                        )

                                elif (
                                    max(s0[1], s1[1]) < y0_in
                                    and max(r0[0], r1[0]) < x1_in
                                ):
                                    return True, "down", rect_in_id

                                elif (
                                    min(s0[1], s1[1]) > y0_in
                                    and min(r0[0], r1[0]) > x1_in
                                ):
                                    return True, "right", rect_in_id
                                else:
                                    return False, "no proper corner intersection", None
                            else:
                                return False, "no proper corner intersection", None

                        elif max_x > x1_in and max_y > y1_in:
                            # --------------UP RIGHT-------------
                            bool_s, intersection_s = intersection_hole_edge(
                                hole, (x1_in, "x"), self.method
                            )
                            bool_r, intersection_r = intersection_hole_edge(
                                hole, (y1_in, "y"), self.method
                            )
                            if bool_r and bool_s:
                                s0, s1 = intersection_s
                                r0, r1 = intersection_r
                                if s0[1] < y1_in < s1[1] or s0[1] > y1_in > s1[1]:
                                    if r0[0] < x1_in < r1[0] or r0[0] > x1_in > r1[0]:
                                        return True, "up_right", rect_in_id
                                    else:
                                        return (
                                            False,
                                            "no proper corner intersection",
                                            None,
                                        )

                                elif (
                                    max(r0[0], r1[0]) < x1_in
                                    and min(s0[1], s1[1]) > y1_in
                                ):
                                    return True, "up", rect_in_id

                                elif (
                                    min(r0[0], r1[0]) > x1_in
                                    and max(s0[1], s1[1]) < y1_in
                                ):
                                    return True, "right", rect_in_id
                                else:
                                    return False, "no proper corner intersection", None

                            else:
                                return False, "no proper corner intersection", None

                        else:
                            return False, "Doesnt intersects with a corner", None
                    else:
                        return (
                            False,
                            "Intersection with too many edges, hole cannot be added",
                            None,
                        )
                else:
                    # ----------------NO INTERSECTIONS --------------------
                    # The hole is completely inside of rect_in
                    return True, None, rect_in_id

    def __add_hole_on_edge(self, intersection_loc, loc: int, hole, convex_hull_hole):
        """
        Args:
            loc:
                Index of the intersected rectangle.
            
        """
        # TODO make it possible to choose how to compute intersection points
        dhole = hole[0]
        x0_in = self.rects_in[loc]["x0"]
        y0_in = self.rects_in[loc]["y0"]
        x1_in = self.rects_in[loc]["x1"]
        y1_in = self.rects_in[loc]["y1"]

        x0_out = self.rects_in[0]["x0"]
        y0_out = self.rects_in[0]["y0"]
        x1_out = self.rects_in[0]["x1"]
        y1_out = self.rects_in[0]["y1"]

        convex_hull_reflected_hole = None
        if self.periodic_boundary:
            convex_hull_reflected_hole = (
                None
                if loc > 0
                else self.__reflect_hole(convex_hull_hole, intersection_loc)
            )
            if convex_hull_reflected_hole is not None:
                if not self.__valid_input(convex_hull_reflected_hole):
                    return False

        # define the intersecting edge
        if intersection_loc == "left":
            edge = (x0_in, "x") if loc > 0 else (x0_out, "x")
        elif intersection_loc == "up":
            edge = (y1_in, "y") if loc > 0 else (y1_out, "y")
        elif intersection_loc == "right":
            edge = (x1_in, "x") if loc > 0 else (x1_out, "x")
        else:
            edge = (y0_in, "y") if loc > 0 else (y0_out, "y")

        bool__, intersections = intersection_hole_edge(hole, edge, self.method)

        if bool__:
            s0, s1 = intersections
            if intersection_loc == "left":
                s0, s1 = (s0, s1) if s0[1] < s1[1] else (s1, s0)

                if loc > 0:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[0] > x0_in
                    )
                    self.rects_in[loc]["edgeLeft"].append([[s0, s1], hole_in, hole_out])

                else:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[0] > x0_out, [0]
                    )
                    if round(hole_in[-1][0], 9) == x0_out:
                        hole_in = hole_in[:-1]
                    elif round(hole_in[0][0], 9) == x0_out:
                        hole_in = hole_in[1:]
                    self.edgeLeft_out.append([s0] + hole_in + [s1])
                    if self.periodic_boundary:
                        self.edgeRight_out.append(
                            [(s1[0] + self.dx, s1[1])]
                            + hole_out
                            + [(s0[0] + self.dx, s0[1])]
                        )

            elif intersection_loc == "up":
                s0, s1 = (s0, s1) if s0[0] < s1[0] else (s1, s0)

                if loc > 0:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[1] < y1_in
                    )
                    self.rects_in[loc]["edgeUp"].append([[s0, s1], hole_in, hole_out])
                else:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[1] < y1_out, [1]
                    )
                    if round(hole_in[-1][1], 9) == y1_out:
                        hole_in = hole_in[:-1]
                    elif round(hole_in[0][1], 9) == y1_out:
                        hole_in = hole_in[1:]

                    self.edgeUp_out.append([s0] + hole_in + [s1])
                    if self.periodic_boundary:
                        self.edgeDown_out.append(
                            [(s1[0], s1[1] - self.dy)]
                            + hole_out
                            + [(s0[0], s0[1] - self.dy)]
                        )

            elif intersection_loc == "right":
                s0, s1 = (s0, s1) if s0[1] > s1[1] else (s1, s0)

                if loc > 0:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[0] < x1_in
                    )
                    self.rects_in[loc]["edgeRight"].append(
                        [[s0, s1], hole_in, hole_out]
                    )
                else:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[0] < x1_out, [2]
                    )
                    if round(hole_in[-1][0], 9) == x1_out:
                        hole_in = hole_in[:-1]

                    elif round(hole_in[0][0], 9) == x1_out:
                        hole_in = hole_in[1:]

                    self.edgeRight_out.append([s0] + hole_in + [s1])
                    if self.periodic_boundary:
                        self.edgeLeft_out.append(
                            [(s1[0] - self.dx, s1[1])]
                            + hole_out
                            + [(s0[0] - self.dx, s0[1])]
                        )

            else:
                s0, s1 = (s0, s1) if s0[0] > s1[0] else (s1, s0)

                if loc > 0:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[1] > y0_in
                    )
                    self.rects_in[loc]["edgeDown"].append([[s0, s1], hole_in, hole_out])
                else:
                    hole_in, hole_out = self.__process_hole(
                        hole, lambda x: x[1] > y0_out, [3]
                    )
                    if round(hole_in[-1][1], 9) == self.y0_out:
                        hole_in = hole_in[:-1]

                    elif round(hole_in[0][1], 9) == self.y0_out:
                        hole_in = hole_in[1:]

                    self.edgeDown_out.append([s0] + hole_in + [s1])
                    if self.periodic_boundary:
                        self.edgeUp_out.append(
                            [(s1[0], s1[1] + self.dy)]
                            + hole_out
                            + [(s0[0], s0[1] + self.dy)]
                        )

            if convex_hull_reflected_hole is not None:
                self.all_holes_conv_hull.append(convex_hull_reflected_hole)
        else:
            return False

    def __add_hole_on_corner(self, intersection_loc, loc, hole):
        # TODO make it possible to choose how to compute intersection points
        dhole = hole[0]
        x0_in = self.rects_in[loc]["x0"]
        y0_in = self.rects_in[loc]["y0"]
        x1_in = self.rects_in[loc]["x1"]
        y1_in = self.rects_in[loc]["y1"]

        x0_out = self.x0_out
        y0_out = self.y0_out
        x1_out = self.x1_out
        y1_out = self.y1_out

        if intersection_loc == "down_left":
            edge_0 = (y0_in, "y") if loc > 0 else (y0_out, "y")
            edge_1 = (x0_in, "x") if loc > 0 else (x0_out, "x")

            bool__, intersections = intersection_hole_corner(
                hole, edge_0, edge_1, intersection_loc, self.method
            )
            s0, s1 = intersections
            if bool__:
                if loc > 0:
                    self.__add_hole_on_corner_in(
                        hole, intersections, intersection_loc, loc
                    )
                else:
                    orientation = lambda x: x[0] > self.x0_out and x[1] > self.y0_out
                    hole_d_l, hole_u_l, hole_u_r, hole_d_r = self.__process_hole(
                        hole, orientation, [3, 0]
                    )

                    if self.periodic_boundary:
                        if not self.__valid_reflection_corner(
                            hole_d_l, hole_u_l, hole_u_r, hole_d_r
                        ):
                            print(
                                "No valid reflection corner, possible solution: more discretization points"
                            )
                            return False

                        s0_x0_y0 = s0
                        s1_x0_y0 = s1

                        s0_x1_y0 = (self.x1_out, s1[1])
                        s1_x1_y0 = ((hole_u_l[-1][0] + hole_d_l[0][0]) / 2, self.y0_out)

                        s0_x0_y1 = (self.x0_out, (hole_d_l[-1][1] + hole_d_r[0][1]) / 2)
                        s1_x0_y1 = (s0[0], self.y1_out)

                        s0_x1_y1 = (s1_x1_y0[0], self.y1_out)
                        s1_x1_y1 = (self.x1_out, s0_x0_y1[1])

                        self.edgeLeft_out.pop(0)
                        self.edgeUp_out.pop(0)
                        self.edgeRight_out.pop(0)
                        self.edgeDown_out.pop(0)

                        self.edgeLeft_out.append([s0_x0_y0] + hole_u_r + [s1_x0_y0])
                        self.edgeUp_out.append([s0_x0_y1] + hole_d_r + [s1_x0_y1])
                        self.edgeRight_out.append([s0_x1_y1] + hole_d_l + [s1_x1_y1])
                        self.edgeDown_out.append([s0_x1_y0] + hole_u_l + [s1_x1_y0])
                    else:
                        self.edgeDown_out.pop(0)
                        self.edgeLeft_out.append([s0] + hole_u_r + [s1])
            else:
                return False

        elif intersection_loc == "left_up":
            edge_0 = (x0_in, "x") if loc > 0 else (x0_out, "x")
            edge_1 = (y1_in, "y") if loc > 0 else (y1_out, "y")
            bool__, intersections = intersection_hole_corner(
                hole, edge_0, edge_1, intersection_loc, self.method
            )
            s0, s1 = intersections
            if bool__:
                if loc > 0:
                    self.__add_hole_on_corner_in(
                        hole, intersections, intersection_loc, loc
                    )
                else:
                    orientation = lambda x: x[0] > self.x0_out and x[1] < self.y1_out
                    hole_d_l, hole_u_l, hole_u_r, hole_d_r = self.__process_hole(
                        hole, orientation, [0, 1]
                    )

                    if self.periodic_boundary:
                        if not self.__valid_reflection_corner(
                            hole_d_l, hole_u_l, hole_u_r, hole_d_r
                        ):
                            return False

                        s0_x1_y1 = (
                            self.x1_out
                            - abs(
                                (
                                    hole_d_l[0][0]
                                    - self.x1_out
                                    + hole_u_l[-1][0]
                                    - self.x1_out
                                )
                                / 2
                            ),
                            self.y1_out,
                        )

                        s1_x1_y1 = (self.x1_out, s0[1])
                        s0_x0_y0 = (s1[0], self.y0_out)
                        s1_x0_y0 = (
                            self.x0_out,
                            abs(hole_u_r[-1][1] + hole_u_l[0][1]) / 2,
                        )
                        s0_x1_y0 = (self.x1_out, s1_x0_y0[1])
                        s1_x1_y0 = (s0_x1_y1[0], self.y0_out)

                        self.edgeUp_out.pop(0)
                        self.edgeRight_out.pop(0)
                        self.edgeLeft_out.pop(0)
                        self.edgeDown_out.pop(0)

                        self.edgeLeft_out.append([s0_x0_y0] + hole_u_r + [s1_x0_y0])
                        self.edgeUp_out.append([s0] + hole_d_r + [s1])
                        self.edgeDown_out.append([s0_x1_y0] + hole_u_l + [s1_x1_y0])
                        self.edgeRight_out.append([s0_x1_y1] + hole_d_l + [s1_x1_y1])
                    else:
                        self.edgeLeft_out.pop(0)
                        self.edgeUp_out.append([s0] + hole_d_r + [s1])
            else:
                return False
        elif intersection_loc == "up_right":
            edge_0 = (y1_in, "y") if loc > 0 else (y1_out, "y")
            edge_1 = (x1_in, "x") if loc > 0 else (x1_out, "x")

            bool__, intersections = intersection_hole_corner(
                hole, edge_0, edge_1, intersection_loc, self.method
            )
            s0, s1 = intersections
            if bool__:
                if loc > 0:
                    self.__add_hole_on_corner_in(
                        hole, intersections, intersection_loc, loc
                    )
                else:
                    orientation = lambda x: x[0] < self.x1_out and x[1] < self.y1_out
                    hole_d_l, hole_u_l, hole_u_r, hole_d_r = self.__process_hole(
                        hole, orientation, [1, 2]
                    )

                    if self.periodic_boundary:
                        if not self.__valid_reflection_corner(
                            hole_d_l, hole_u_l, hole_u_r, hole_d_r
                        ):
                            return False

                        s0_x1_y1 = s0
                        s1_x1_y1 = s1

                        s0_x1_y0 = (self.x1_out, (hole_u_l[0][1] + hole_u_r[-1][1]) / 2)
                        s1_x1_y0 = (s0[0], self.y0_out)

                        s0_x0_y1 = (self.x0_out, s1[1])
                        s1_x0_y1 = ((hole_u_r[0][0] + hole_d_r[-1][0]) / 2, self.y1_out)

                        s0_x0_y0 = (s1_x0_y1[0], self.y0_out)
                        s1_x0_y0 = (self.x0_out, s0_x1_y0[1])

                        self.edgeLeft_out.pop(0)
                        self.edgeUp_out.pop(0)
                        self.edgeRight_out.pop(0)
                        self.edgeDown_out.pop(0)

                        self.edgeLeft_out.append([s0_x0_y0] + hole_u_r + [s1_x0_y0])
                        self.edgeUp_out.append([s0_x0_y1] + hole_d_r + [s1_x0_y1])
                        self.edgeDown_out.append([s0_x1_y0] + hole_u_l + [s1_x1_y0])
                        self.edgeRight_out.append([s0_x1_y1] + hole_d_l + [s1_x1_y1])
                    else:
                        self.edgeUp_out.pop(0)
                        self.edgeRight_out.append([s0] + hole_d_l + [s1])

            else:
                return False
        else:
            edge_0 = (x1_in, "x") if loc > 0 else (x1_out, "x")
            edge_1 = (y0_in, "y") if loc > 0 else (y0_out, "y")

            bool__, intersections = intersection_hole_corner(
                hole, edge_0, edge_1, intersection_loc, self.method
            )
            s0, s1 = intersections
            if bool__:
                if loc > 0:
                    self.__add_hole_on_corner_in(
                        hole, intersections, intersection_loc, loc
                    )
                else:
                    orientation = lambda x: x[0] < self.x1_out and x[1] > self.y0_out
                    hole_d_l, hole_u_l, hole_u_r, hole_d_r = self.__process_hole(
                        hole, orientation, [2, 3]
                    )

                    if self.periodic_boundary:
                        if not self.__valid_reflection_corner(
                            hole_d_l, hole_u_l, hole_u_r, hole_d_r
                        ):
                            return False

                        s0, s1 = intersections

                        s0_x1_y0 = s0
                        s1_x1_y0 = s1

                        s0_x0_y1 = (self.x0_out, (hole_d_r[0][1] + hole_d_l[-1][1]) / 2)
                        s1_x0_y1 = ((hole_u_r[0][0] + hole_d_r[-1][0]) / 2, self.y1_out)

                        s0_x1_y1 = (s1_x1_y0[0], self.y1_out)
                        s1_x1_y1 = (self.x1_out, s0_x0_y1[1])

                        s0_x0_y0 = (s1_x0_y1[0], self.y0_out)
                        s1_x0_y0 = (self.x0_out, s0_x1_y0[1])

                        self.edgeLeft_out.pop(0)
                        self.edgeUp_out.pop(0)
                        self.edgeRight_out.pop(0)
                        self.edgeDown_out.pop(0)

                        self.edgeLeft_out.append([s0_x0_y0] + hole_u_r + [s1_x0_y0])
                        self.edgeUp_out.append([s0_x0_y1] + hole_d_r + [s1_x0_y1])
                        self.edgeDown_out.append([s0_x1_y0] + hole_u_l + [s1_x1_y0])
                        self.edgeRight_out.append([s0_x1_y1] + hole_d_l + [s1_x1_y1])
                    else:
                        self.edgeRight_out.pop(0)
                        self.edgeDown_out.append([s0] + hole_u_l + [s1])

            else:
                return False

    def __add_hole_on_corner_in(self, hole, intersections, position, rect_id):
        x0_in = self.rects_in[rect_id]["x0"]
        y0_in = self.rects_in[rect_id]["y0"]
        x1_in = self.rects_in[rect_id]["x1"]
        y1_in = self.rects_in[rect_id]["y1"]

        x0_out = self.x0_out
        y0_out = self.y0_out
        x1_out = self.x1_out
        y1_out = self.y1_out
        s0, s1 = intersections
        if position == "down_left":
            orientation_corner_case = (lambda x: x[0] < x0_in, lambda x: x[1] > y0_in)
            hole_in, hole_out = self.__process_hole(
                hole,
                lambda x: x[0] > x0_in and x[1] > y0_in,
                orientation_corner_case=orientation_corner_case,
            )
            self.rects_in[rect_id]["edgeDown"].pop(0)
            self.rects_in[rect_id]["edgeLeft"].append([[s0, s1], hole_in, hole_out])
        elif position == "left_up":
            orientation_corner_case = (lambda x: x[0] > x0_in, lambda x: x[1] > y1_in)
            hole_in, hole_out = self.__process_hole(
                hole,
                lambda x: x[0] > x0_in and x[1] < y1_in,
                orientation_corner_case=orientation_corner_case,
            )
            self.rects_in[rect_id]["edgeLeft"].pop(0)
            self.rects_in[rect_id]["edgeUp"].append([[s0, s1], hole_in, hole_out])
        elif position == "up_right":
            orientation_corner_case = (lambda x: x[0] > x1_in, lambda x: x[1] < y1_in)
            hole_in, hole_out = self.__process_hole(
                hole,
                lambda x: x[0] < x1_in and x[1] < y1_in,
                orientation_corner_case=orientation_corner_case,
            )
            self.rects_in[rect_id]["edgeUp"].pop(0)
            self.rects_in[rect_id]["edgeRight"].append([[s0, s1], hole_in, hole_out])

        else:
            orientation_corner_case = (lambda x: x[1] < y0_in, lambda x: x[0] > x1_in)
            hole_in, hole_out = self.__process_hole(
                hole,
                lambda x: x[0] < x1_in and x[1] > y0_in,
                orientation_corner_case=orientation_corner_case,
            )
            self.rects_in[rect_id]["edgeRight"].pop(0)
            self.rects_in[rect_id]["edgeDown"].append([[s0, s1], hole_in, hole_out])

    def edges_out(self):
        return [
            sorted(self.edgeLeft_out, key=lambda x: x[0][1]),
            sorted(self.edgeUp_out, key=lambda x: x[0][0]),
            sorted(self.edgeRight_out, key=lambda x: -x[0][1]),
            sorted(self.edgeDown_out, key=lambda x: -x[0][0]),
        ]

    def edges_in(self, rect_id):
        """
        sorts the holes on each edge of rect_id.
        sort the holes clockwise
        for x in self.rects_in[rect_id]["someEdge]:
            x[0] = intersection points (can be of length 1 or 2, 1 if hole is vertex)
            x[0][0] = first intersection point = tuple (x-coordinate, y-coordinate)
        """
        return [
            sorted(
                self.rects_in[rect_id]["edgeLeft"],
                key=lambda x: x[0][0][1] if len(x) == 3 else x[0][1],
            ),
            sorted(
                self.rects_in[rect_id]["edgeUp"],
                key=lambda x: x[0][0][0] if len(x) == 3 else x[0][0],
            ),
            sorted(
                self.rects_in[rect_id]["edgeRight"],
                key=lambda x: -x[0][0][1] if len(x) == 3 else -x[0][1],
            ),
            sorted(
                self.rects_in[rect_id]["edgeDown"],
                key=lambda x: -x[0][0][0] if len(x) == 3 else -x[0][0],
            ),
        ]

    def __process_hole(
        self, hole, orientation, ref_dir=None, orientation_corner_case=None
    ):
        # orientation(hole[i]) = True if in else False w.r.t. an edge or an corner
        # hole = [p0, p1, ..., pN]
        # ref_dir tells where to shift the hole in the periodic case
        # TODO: Check this
        dhole = hole[0]
        bool_ = False
        counter = 0
        maxCounter = len(dhole)
        # shift the hole as long as the last point is in and the first one is out
        # possible problem:
        #   -hole intersects only with boundary!

        # TODO hole  must be discretization and parametric information aswell
        #
        #
        #    2 Faelle :
        #
        #  (1) hole is nach verfeinerung gar nicht im schnitt mit dem rechteck  :
        #      # gib information  außen / innen zurück sodass callende methode die richtige methode aufrufen kann
        #      # zum beispiel optimierung oder  N punkte zwischen den "außen punkten" nehmen und testen ob einer im rechteck liegt. N >> 1
        #
        #
        #  (2) hole schneidet das gebiet :
        #      mache den while loop unten  öfter :
        #
        #      2.1  hole ist kein stella : abbruch nach einem loop mit could not proccess hole, rise Notimplemented Expection ("cannot generate point in between")
        #      2.2  hole ist ein stella :
        #               for dem ersten shift, speicher das winkel array mit was zu den tupeln gehört und shifte es immer mit
        #
        #           kriege raus welche 2 punkte den linearen schnitt mit dem Rechteck verursachen,  neheme die dazugehörigen winkel ( modulu 2pi)
        #           und füge nen neuen winkel in der mitte in das winkel array und das tupel array
        #
        #           Falls der hinzugefügte punkt im inneren liegt und der orientation test also geht fertig
        #           sonst füge weitere punkte hinzu .
        #
        #
        #
        #      füge nen punkt hinzu
        #           starte while loop  bis es funktioniert
        #
        #
        #
        #

        while not bool_ and counter <= maxCounter + 2:
            dhole = shift(dhole, 1)
            bool_ = orientation(dhole[-1]) and not orientation(dhole[0])
            counter += 1

        if not bool_:
            # In this case the polygon dhole intersects the rect, but none of the discretization points is inside the rect.
            # This means all dpoint are outside => make inner loop only with intersection points and corner point
            #   Problem: still need to shift outer loop
            shift_index = None
            for i in range(maxCounter):
                p = dhole[i]
                p_prev = dhole[i - 1] if i > 0 else dhole[maxCounter - 1]
                if orientation_corner_case[0](p) and orientation_corner_case[1](p):
                    if (not orientation_corner_case[0](p_prev)) and (
                        not orientation_corner_case[1](p_prev)
                    ):
                        shift_index = i
            if shift_index is None:
                print("Should not enter this case.")
                exit()

            return [], shift(dhole, shift_index)

        hole_in = [point for point in dhole if orientation(point)]

        if ref_dir is None:
            hole_out = [point for point in dhole if not orientation(point)]
            return hole_in, hole_out

        elif len(ref_dir) == 1:
            # edge case
            if ref_dir[0] == 2:
                hole_out = [(p[0] - self.dx, p[1]) for p in dhole if not orientation(p)]

            elif ref_dir[0] == 3:
                hole_out = [(p[0], p[1] + self.dy) for p in dhole if not orientation(p)]

            elif ref_dir[0] == 0:
                hole_out = [(p[0] + self.dx, p[1]) for p in dhole if not orientation(p)]

            else:
                hole_out = [(p[0], p[1] - self.dy) for p in dhole if not orientation(p)]

            return hole_in, hole_out

        elif len(ref_dir) == 2:
            # corner case
            if ref_dir[0] == 0 and ref_dir[1] == 1:
                # up/left
                hole_d_r = hole_in
                hole_d_l = [
                    (p[0] + self.dx, p[1])
                    for p in dhole
                    if p[0] < self.x0_out and p[1] < self.y1_out
                ]
                hole_u_l = [
                    (p[0] + self.dx, p[1] - self.dy)
                    for p in dhole
                    if p[0] < self.x0_out and p[1] > self.y1_out
                ]
                hole_u_r = [
                    (p[0], p[1] - self.dy)
                    for p in dhole
                    if p[0] > self.x0_out and p[1] > self.y1_out
                ]

            elif ref_dir[0] == 1 and ref_dir[1] == 2:
                # up/right
                hole_u_r = [
                    (p[0] - self.dx, p[1] - self.dy)
                    for p in dhole
                    if p[0] > self.x1_out and p[1] > self.y1_out
                ]
                hole_u_l = [
                    (p[0], p[1] - self.dy)
                    for p in dhole
                    if p[0] < self.x1_out and p[1] > self.y1_out
                ]

                hole_d_r = [
                    (p[0] - self.dx, p[1])
                    for p in dhole
                    if p[0] > self.x1_out and p[1] < self.y1_out
                ]
                hole_d_l = hole_in

            elif ref_dir[0] == 2 and ref_dir[1] == 3:
                # down/right
                hole_u_r = [
                    (p[0] - self.dx, p[1])
                    for p in dhole
                    if p[0] > self.x1_out and p[1] > self.y0_out
                ]
                hole_u_l = hole_in
                hole_d_r = [
                    (p[0] - self.dx, p[1] + self.dy)
                    for p in dhole
                    if p[0] > self.x1_out and p[1] < self.y0_out
                ]
                hole_d_l = [
                    (p[0], p[1] + self.dy)
                    for p in dhole
                    if p[0] < self.x1_out and p[1] < self.y0_out
                ]

            else:
                # down/left
                hole_u_r = hole_in
                hole_u_l = [
                    (p[0] + self.dx, p[1])
                    for p in dhole
                    if p[0] < self.x0_out and p[1] > self.y0_out
                ]
                hole_d_r = [
                    (p[0], p[1] + self.dy)
                    for p in dhole
                    if p[0] > self.x0_out and p[1] < self.y0_out
                ]
                hole_d_l = [
                    (p[0] + self.dx, p[1] + self.dy)
                    for p in dhole
                    if p[0] < self.x0_out and p[1] < self.y0_out
                ]

            return hole_d_l, hole_u_l, hole_u_r, hole_d_r

    def __reflect_hole(self, convex_hull_hole, intersection_loc):
        if intersection_loc == "left":
            return [(x[0] + self.dx, x[1]) for x in convex_hull_hole]
        elif intersection_loc == "up":
            return [(x[0], x[1] - self.dy) for x in convex_hull_hole]
        elif intersection_loc == "right":
            return [(x[0] - self.dx, x[1]) for x in convex_hull_hole]
        else:
            return [(x[0], x[1] + self.dy) for x in convex_hull_hole]

    def get_geometry(self, lc, lc_subrects, filled=False):
        """
        I)  Write boundary of rect_0 (Main rect that contains all other rects)
        II) Write boundaries of rect_1,...,rect_n
            - for each rect write inner and outer boundary
                - rect_i has inner loop ID 2*i and outer loop ID 2*i + 1
                - outer loop is needed for surface of rect_0
                - inner loop is needed for surface of rect_i
        III) Write all line loops of holes
            - store line loop ID and the rect containing it, e.g. (rect_ID, loop_ID)
        VI) If holes need to be filled:
            - write line loops for holes on boundaries
            - for every hole loop make surface

        """

        geometry = {
            "points": {},
            "lines": {},
            "lineLoops": {},
            "surfaces": {},
            "physicalSurfaces": {Id: [Id] for Id in range(1, self.nRects + 2)},
        }
        freePointID = 1
        freeLineID = 1
        freeLineLoopID = (
            2 * (len(self.rects_in) - 1) + 2
        )  # lower IDs are reserved for rect loops
        freeSurfaceID = len(self.rects_in) + 1

        # ===================================
        #   I) Write boundary of rect_0 (Main rect that contains all other rects)
        # ===================================
        boundaryLoopsOut = []
        edges_out = self.edges_out()
        firstLineID = freeLineID
        firstPointID = freePointID
        lastPointID = None
        for i in range(4):
            edge = edges_out[i]
            nholes = len(edge)
            for j in range(nholes):
                hole = edge[j]
                if len(hole) == 1:
                    # vertex case
                    x, y = hole[0]
                    coords = (x, y, 0)
                    geometry["points"][freePointID] = {
                        "id": freePointID,
                        "coords": coords,
                        "lc": lc,
                    }
                    freePointID += 1

                    if lastPointID is not None:
                        geometry["lines"][freeLineID] = (lastPointID, freePointID - 1)
                        freeLineID += 1
                    lastPointID = freePointID - 1

                else:
                    # hole case
                    nPoints = len(hole)
                    lPoints = []
                    for k in range(nPoints):
                        x, y = hole[k]
                        coords = (x, y, 0)
                        geometry["points"][freePointID] = {
                            "id": freePointID,
                            "coords": coords,
                            "lc": lc,
                        }
                        freePointID += 1
                        lPoints.append(freePointID - 1)

                    firstIntersectionPoint = freePointID - nPoints

                    if lastPointID is not None:
                        geometry["lines"][freeLineID] = (
                            lastPointID,
                            firstIntersectionPoint,
                        )
                        freeLineID += 1
                    lastPointID = freePointID - 1
                    # make lines
                    loop = []
                    nLines = nPoints - 1
                    for k in range(nLines):
                        start = lPoints[k]
                        end = lPoints[k + 1]
                        geometry["lines"][freeLineID] = (start, end)
                        freeLineID += 1
                        loop.append(freeLineID - 1)
                    boundaryLoopsOut.append(loop)

        geometry["lines"][freeLineID] = (lastPointID, firstPointID)
        freeLineID += 1
        geometry["lineLoops"][1] = list(range(firstLineID, freeLineID))

        # ===================================
        #   II) Write boundaries of rect_1,...,rect_n
        # ===================================
        boundaryLoops = []
        for rect_ID, rect_in in self.rects_in.items():
            if rect_ID > 0:
                inner_boundary = []
                outer_boundary = []
                edges_in = self.edges_in(rect_ID)
                firstPointID = None
                lastPointID = None
                for i in range(4):
                    edge = edges_in[i]
                    nholes = len(edge)
                    for j in range(nholes):
                        if i == 0 and j == 0:
                            firstPointID = freePointID
                        hole = edge[j]
                        if len(hole) == 1:
                            # vertex case
                            coords = (hole[0][0], hole[0][1], 0)
                            geometry["points"][freePointID] = {
                                "id": freePointID,
                                "coords": coords,
                                "lc": lc_subrects,
                            }
                            freePointID += 1
                            if lastPointID is not None:
                                geometry["lines"][freeLineID] = (
                                    lastPointID,
                                    freePointID - 1,
                                )
                                freeLineID += 1
                                inner_boundary.append(freeLineID - 1)
                                outer_boundary.append(freeLineID - 1)

                            lastPointID = freePointID - 1
                        else:
                            intersections = hole[0]
                            hole_in = hole[1]
                            hole_out = hole[2]

                            firstIntersectionPoint = freePointID
                            secondIntersectionPoint = freePointID + 1

                            if i == 1 and j == 0 and len(edges_in[0]) == 0:
                                firstPointID = firstIntersectionPoint

                            # write intersection points
                            coords = (intersections[0][0], intersections[0][1], 0)
                            geometry["points"][freePointID] = {
                                "id": freePointID,
                                "coords": coords,
                                "lc": lc_subrects,
                            }
                            freePointID += 1

                            coords = (intersections[1][0], intersections[1][1], 0)
                            geometry["points"][freePointID] = {
                                "id": freePointID,
                                "coords": coords,
                                "lc": lc_subrects,
                            }
                            freePointID += 1

                            # write points of inner hole
                            nPointsIn = len(hole_in)
                            lPointsIn = []
                            for k in range(nPointsIn):
                                point = hole_in[k]
                                x, y = point
                                coords = (x, y, 0)
                                geometry["points"][freePointID] = {
                                    "id": freePointID,
                                    "coords": coords,
                                    "lc": lc_subrects,
                                }
                                freePointID += 1
                                lPointsIn.append(freePointID - 1)
                            lPointsIn.append(secondIntersectionPoint)
                            lPointsIn = [firstIntersectionPoint] + lPointsIn
                            # write points of outer hole
                            nPointsOut = len(hole_out)
                            lPointsOut = []
                            for k in range(nPointsOut):
                                point = hole_out[k]
                                x, y = point
                                coords = (x, y, 0)
                                geometry["points"][freePointID] = {
                                    "id": freePointID,
                                    "coords": coords,
                                    "lc": lc_subrects,
                                }
                                freePointID += 1
                                lPointsOut.append(freePointID - 1)
                            lPointsOut.reverse()
                            lPointsOut.append(secondIntersectionPoint)
                            lPointsOut = [firstIntersectionPoint] + lPointsOut
                            # write connection lines
                            if lastPointID is not None:
                                geometry["lines"][freeLineID] = (
                                    lastPointID,
                                    firstIntersectionPoint,
                                )
                                freeLineID += 1
                                inner_boundary.append(freeLineID - 1)
                                outer_boundary.append(freeLineID - 1)

                            lastPointID = secondIntersectionPoint

                            # write lines corresponding to inner points
                            linesIn = []
                            nLinesIn = nPointsIn + 1
                            for k in range(nLinesIn):
                                geometry["lines"][freeLineID] = (
                                    lPointsIn[k],
                                    lPointsIn[k + 1],
                                )
                                freeLineID += 1
                                inner_boundary.append(freeLineID - 1)
                                linesIn.append(freeLineID - 1)

                            # write lines corresponding to outer points
                            linesOut = []
                            nLinesOut = nPointsOut + 1
                            for k in range(nLinesOut):
                                geometry["lines"][freeLineID] = (
                                    lPointsOut[k],
                                    lPointsOut[k + 1],
                                )
                                freeLineID += 1
                                outer_boundary.append(freeLineID - 1)
                                linesOut.append(freeLineID - 1)

                            boundaryLoops.append((rect_ID, linesIn, linesOut))

                geometry["lines"][freeLineID] = (lastPointID, firstPointID)
                freeLineID += 1
                inner_boundary.append(freeLineID - 1)
                outer_boundary.append(freeLineID - 1)

                geometry["lineLoops"][2 * rect_ID] = inner_boundary
                geometry["lineLoops"][2 * rect_ID + 1] = outer_boundary

        # ===================================
        #   III)
        # ===================================
        lInnerHoles = [[] for i in self.rects_in]
        for hole in self.holes_between:
            nPoints = len(hole)
            lPointID = []
            lineLoop = []
            for point in hole:
                x, y = point
                coords = (x, y, 0)
                geometry["points"][freePointID] = {
                    "id": freePointID,
                    "coords": coords,
                    "lc": None,
                }
                freePointID += 1
                lPointID.append(freePointID - 1)

            nPoints = len(hole)
            for i in range(nPoints):
                start = lPointID[i]
                end = lPointID[i + 1] if i < nPoints - 1 else lPointID[0]
                geometry["lines"][freeLineID] = (start, end)
                freeLineID += 1
                lineLoop.append(freeLineID - 1)

            geometry["lineLoops"][freeLineLoopID] = lineLoop
            freeLineLoopID += 1
            lInnerHoles[0].append(freeLineLoopID - 1)

        for rect_ID, rect_in in self.rects_in.items():
            if rect_ID > 0:
                holes = rect_in["holes_in"]
                for hole in holes:
                    nPoints = len(hole)
                    lPointID = []
                    lineLoop = []
                    for point in hole:
                        x, y = point
                        coords = (x, y, 0)
                        geometry["points"][freePointID] = {
                            "id": freePointID,
                            "coords": coords,
                            "lc": None,
                        }
                        freePointID += 1
                        lPointID.append(freePointID - 1)

                    nPoints = len(hole)
                    for i in range(nPoints):
                        start = lPointID[i]
                        end = lPointID[i + 1] if i < nPoints - 1 else lPointID[0]
                        geometry["lines"][freeLineID] = (start, end)
                        freeLineID += 1
                        lineLoop.append(freeLineID - 1)

                    geometry["lineLoops"][freeLineLoopID] = lineLoop
                    freeLineLoopID += 1
                    lInnerHoles[rect_ID].append(freeLineLoopID - 1)

        # ==================================
        #   IV)
        # ==================================
        boundaryHolesSurfaces = [[] for i in self.rects_in]
        lboundaryHoles = [[] for i in self.rects_in]
        if filled:
            for loop in boundaryLoops:
                rect_ID, loopIn, loopOut = loop
                point0_ID = geometry["lines"][loopIn[0]][0]
                point1_ID = geometry["lines"][loopIn[-1]][1]

                x0, y0, z0 = geometry["points"][point0_ID]["coords"]
                x1, y1, z1 = geometry["points"][point1_ID]["coords"]

                if x0 - x1 != 0 and y0 - y1 != 0:
                    # corner case, we need to add a point
                    x, y = None, None
                    if x0 == self.rects_in[rect_ID]["x0"]:
                        x = x0
                        y = self.rects_in[rect_ID]["y1"]

                    elif y0 == self.rects_in[rect_ID]["y1"]:
                        # up right
                        x = self.rects_in[rect_ID]["x1"]
                        y = y0
                    elif x0 == self.rects_in[rect_ID]["x1"]:
                        # down right
                        x = x0
                        y = self.rects_in[rect_ID]["y0"]
                    elif y0 == self.rects_in[rect_ID]["y0"]:
                        # down left
                        x = self.rects_in[rect_ID]["x0"]
                        y = y0
                    else:
                        raise ValueError("Dunno")

                    coords = (x, y, 0)
                    geometry["points"][freePointID] = {
                        "id": freePointID,
                        "coords": coords,
                        "lc": None,
                    }
                    freePointID += 1

                    geometry["lines"][freeLineID] = (point1_ID, freePointID - 1)
                    freeLineID += 1

                    geometry["lines"][freeLineID] = (freePointID - 1, point0_ID)
                    freeLineID += 1

                    geometry["lineLoops"][freeLineLoopID] = loopIn + [
                        freeLineID - 2,
                        freeLineID - 1,
                    ]
                    freeLineLoopID += 1
                    lboundaryHoles[rect_ID].append(freeLineLoopID - 1)

                    geometry["lineLoops"][freeLineLoopID] = loopOut + [
                        freeLineID - 2,
                        freeLineID - 1,
                    ]
                    freeLineLoopID += 1
                    lboundaryHoles[0].append(freeLineLoopID - 1)

                else:
                    geometry["lines"][freeLineID] = (point1_ID, point0_ID)
                    freeLineID += 1

                    geometry["lineLoops"][freeLineLoopID] = loopIn + [freeLineID - 1]
                    freeLineLoopID += 1
                    lboundaryHoles[rect_ID].append(freeLineLoopID - 1)

                    geometry["lineLoops"][freeLineLoopID] = loopOut + [freeLineID - 1]
                    freeLineLoopID += 1
                    lboundaryHoles[0].append(freeLineLoopID - 1)

            for loop in boundaryLoopsOut:
                # first check if we need to add a point

                point0_ID = geometry["lines"][loop[0]][0]
                point1_ID = geometry["lines"][loop[-1]][1]

                x0, y0, z0 = geometry["points"][point0_ID]["coords"]
                x1, y1, z1 = geometry["points"][point1_ID]["coords"]

                if x0 - x1 != 0 and y0 - y1 != 0:
                    # corner case, we need to add a point
                    x, y = None, None
                    if x0 == self.x0_out:
                        x = x0
                        y = self.y1_out

                    elif y0 == self.y1_out:
                        # up right
                        x = self.x1_out
                        y = y0
                    elif x0 == self.x1_out:
                        # down right
                        x = x0
                        y = self.y0_out
                    elif y0 == self.y0_out:
                        # down left
                        x = self.x0_out
                        y = y0
                    else:
                        raise ValueError("Dunno")

                    coords = (x, y, 0)
                    geometry["points"][freePointID] = {
                        "id": freePointID,
                        "coords": coords,
                        "lc": None,
                    }
                    freePointID += 1

                    geometry["lines"][freeLineID] = (point1_ID, freePointID - 1)
                    freeLineID += 1
                    geometry["lines"][freeLineID] = (freePointID - 1, point0_ID)
                    freeLineID += 1

                    geometry["lineLoops"][freeLineLoopID] = loop + [
                        freeLineID - 2,
                        freeLineID - 1,
                    ]
                    freeLineLoopID += 1
                    lboundaryHoles[0].append(freeLineLoopID - 1)

                else:
                    geometry["lines"][freeLineID] = (point1_ID, point0_ID)
                    freeLineID += 1
                    geometry["lineLoops"][freeLineLoopID] = loop + [freeLineID - 1]
                    freeLineLoopID += 1
                    lboundaryHoles[0].append(freeLineLoopID - 1)

        # ==================================
        # Write main surfaces
        # ==================================
        geometry["surfaces"][1] = (
            [1] + [2 * i + 1 for i in range(1, len(self.rects_in))] + lInnerHoles[0]
        )

        for rect_ID in range(1, len(self.rects_in)):
            geometry["surfaces"][rect_ID + 1] = [2 * rect_ID] + lInnerHoles[rect_ID]

        if filled:
            # make surfaces of boundary holes
            lLoops = [
                lInnerHoles[i] + lboundaryHoles[i] for i in range(len(self.rects_in))
            ]

            physicalHoles = [[] for i in range(self.nRects + 1)]
            for Id, loops in enumerate(lLoops):
                for loop in loops:
                    geometry["surfaces"][freeSurfaceID] = [loop]
                    freeSurfaceID += 1
                    physicalHoles[Id].append(freeSurfaceID - 1)

                    # physicalHoles[1].append(freeSurfaceID-1)

            for Id, surfaces in enumerate(physicalHoles):
                tag = Id + 2 + self.nRects
                # print(" ADD PYHSICAL HOLE SURFACES : ")
                # print(" tag = ", str(tag))
                # print("surfaces = ", str(surfaces))
                if len(surfaces) == 0:
                    continue
                geometry["physicalSurfaces"][tag] = surfaces

        return geometry


def shift(array, n):
    # doesnt work with numpy arrays!
    # left shift
    m = len(array)
    n = n % m
    return array[n:] + array[: -(m - n)]
