# TODO: Is this module completely deprecated?
from bubbles.two_d.writeGeo import WriteGeo


class DomainWithHoles(WriteGeo):
    
    def __init__(self, file_name, dir = None):
        raise ValueError("Use of this class is deprecated! Use instead only class topology with method get_geometry()")
        WriteGeo.__init__(self, file_name, dir)
        self.line_loops_in = []
        self.line_loops_out = []
        self.line_loop_inner_rect = [] # Id of inner rect
        self.line_loop_outer_rect = None
        self.physical_in = []
        self.physical_out = []

    def write_files(self, refinements = [(.2, .2)]):
        for refs in refinements:
            lc_in, lc_out = refs
            self.write_file(topology, lc_in, lc_out)

    def write_file(self, topology, lc_in, lc_out, flag):
        assert flag == "filled" or flag == "hole"
        self.write("lc_in = " + str(lc_in))
        self.write("lc_out = " + str(lc_out))
        
        holes_rects = {}
        for rect_ID, rect_in in topology.rects_in.items():
            if rect_ID > 0:
                self.comment("=========================")
                self.comment("Write inner rect:")
                self.comment("=========================")
                # Iterate over inner rects:
                # write points, lines, line loops, surfaces and physical surfaces
                rect_in = topology.rects_in[rect_ID]
                edges_in = topology.edges_in(rect_ID) 
                x0_in = rect_in["x0"]
                y0_in = rect_in["y0"]
                x1_in = rect_in["x1"]
                y1_in = rect_in["y1"]
                rect_points = [(x0_in,y0_in), (x1_in, y1_in)]
                holes_in, holes_out = self.write_boundary_in(edges_in, rect_points)
                rect_loop_id = self.lineLoopId-1
                holes_interior_in = self.write_holes_in(rect_in["holes_in"])
                self.write_plane_surface([rect_loop_id]+holes_in + holes_interior_in)
                self.write_physical_surface(str(rect_ID),str(self.planeSurfaceId-1))

                holes_rects[rect_ID] = holes_in + holes_interior_in

        
        if flag == "filled":
            self.comment("=========================")
            self.comment("Write surfaces of rects")
            self.comment("=========================")
            for rect_ID, holes in holes_rects.items():
                physical_surfaces = []
                for hole in holes:
                    self.write_plane_surface([hole])
                    physical_surfaces.append(str(self.planeSurfaceId-1))
                
                tag = str(rect_ID)+str(rect_ID)
                self.write_physical_surface(tag, physical_surfaces)

        self.comment("=========================")
        self.comment("Write boundary outer rect:")
        self.comment("=========================")

        self.write_boundary_out(topology)

        self.comment("=========================")
        self.comment("Write holes between rects and outer rect:")
        self.comment("=========================")

        self.write_holes_out(topology)

        if flag == "filled":
            for j, hole in enumerate(self.line_loops_out):
                self.write_plane_surface([hole])
                self.physical_out.append(str(self.planeSurfaceId - 1))
               

        self.write_plane_surface([self.line_loop_outer_rect]
         + self.line_loops_out + self.line_loop_inner_rect)
        self.write_physical_surface("-1", str(self.planeSurfaceId - 1))
        self.write_physical_surface("-11", self.physical_out)
       

    def write_holes_in(self, holes_in):
        loops_in = []
        for hole in holes_in:
            self.write_polygon(hole, "lc_in")
            loops_in.append(self.lineLoopId-1)
            self.line_loops_in.append(self.lineLoopId - 1)
        return loops_in

    def write_holes_out(self, topology):
        for hole in topology.holes_between:
            self.write_polygon(hole, "lc_out")#self.lc_out

            if len(self.line_loops_out) == 0:
                self.line_loops_out = [self.lineLoopId - 1] + self.line_loops_out[1:]
            else:
                self.line_loops_out = [self.line_loops_out[0]] + [self.lineLoopId - 1] + self.line_loops_out[1:]

    def write_boundary_in(self, edges_in,rect_points):
        """
        Here we write the inner rectangle, therfore we write two line loops.
        Holes on the boundary are NOT defined "through" it's own line loop.
        """
        for i in range(4):
            print("============")
            print(i)
            print("_________")
            for j in edges_in[i]:
                if len(j) > 2:
                    print(len(j[1]), len(j[2]))
                print(j)
                print("_________")
        print("==============")
        print(edges_in[0][0])

        x0_in,y0_in = rect_points[0]
        x1_in,y1_in = rect_points[1]

        line_loop_in_new = []
        counter = 0
        temp = 0
        first_point_id = self.pointId
        holes_in, holes_out = [],[]
        for i, edge in enumerate(edges_in):
            # We write the boundary in following order:
            # left -> up -> right -> down

            for j, hole in enumerate(edge):
                """
                hole is a list.
                if len(list) == 1 then it is just a vertice of the inner rectangle,
                otherwise it is a hole.
                """
                if len(hole) == 1:
                    self.write_point(hole[0][0], hole[0][1], "lc_in")
                    counter += 1

                    if j > 0 or i > 0:
                        self.write_line(self.pointId - temp - 2, self.pointId - 1)
                        
                        line_loop_in_new.append(self.lineId - 1)

                        temp = 0
                else:
                    s0, s1 = hole[0] # start-  and endpoint of the hole (intersections)
                    hole_in = hole[1]
                    hole_out = hole[2]
                    hole_line_loop_in = []
                    hole_line_loop_out = []
                    # write first intersection point
                    self.write_point(s0[0], s0[1], "lc_in")
                    counter += 1
                    # write connecting line with previous hole

                    if i == 1 and j == 0 and len(edges_in[0]) != 0:
                        self.write_line(self.pointId - 2 - temp, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    elif i > 1 and j == 0:
                        # TODO check bug !
                        self.write_line(self.pointId - 2 - temp, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    elif j > 0:
                        self.write_line(self.pointId - 2 - temp, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    #write points
                    for h in hole_in:
                        self.write_point(h[0], h[1], "lc_in")
                        counter += 1
                    self.write_point(s1[0], s1[1], "lc_in")
                    counter += 1
                    for h in hole_out[1:-1]:
                        self.write_point(h[0], h[1], "lc_in")
                        counter += 1

                    # write lines
                    n_holes_in = len(hole_in)
                    n_holes_out = len(hole_out) - 2
                    n_tot = n_holes_in + n_holes_out + 2
                    temp = n_holes_out # Before + 2

                    # write lines in
                    for k in range(n_holes_in + 1):
                        self.write_line(self.pointId - n_tot + k, self.pointId - n_tot + k + 1)
                        hole_line_loop_in.append(self.lineId - 1)

                    if s0[i % 2] == s1[i % 2]:
                        # case: edge
                        self.write_line(self.pointId - 1 - n_holes_out, self.pointId - n_tot)
                        line_loop_in_new.append(- (self.lineId - 1))
                        hole_line_loop_in.append(self.lineId - 1)
                        hole_line_loop_out.append(-(self.lineId - 1))

                        # write lines out
                        for k in range(n_holes_out):
                            self.write_line(self.pointId - 1 - n_holes_out + k, self.pointId - n_holes_out + k)
                            hole_line_loop_out.append(self.lineId - 1)

                        self.write_line(self.pointId - 1, self.pointId - n_tot)
                        hole_line_loop_out.append(self.lineId - 1)

                    else:
                        # case: corner
                        corner_points = [
                        (x0_in, y0_in),
                        (x0_in, y1_in),
                        (x1_in, y1_in),
                        (x1_in, y0_in),
                        ]
                        p = corner_points[i]
                        self.write_point(p[0], p[1], "lc_in")
                        counter += 1
                        temp += 1
                        # write lines between hole_in and hole_out
                        self.write_line(self.pointId - 1 - n_holes_out -1, self.pointId - 1)

                        hole_line_loop_in.append(self.lineId - 1)

                        self.write_line(self.pointId - 1, self.pointId - 1 - n_tot)
                        line_loop_in_new.append(- (self.lineId - 1))
                        line_loop_in_new.append(- (self.lineId - 2))
                        hole_line_loop_in.append(self.lineId - 1)
                        hole_line_loop_out.append(-(self.lineId - 1))
                        hole_line_loop_out.append(-(self.lineId - 2))

                        # write lines out
                        for k in range(n_holes_out):
                            self.write_line(self.pointId - 2 - n_holes_out + k, self.pointId - 1 - n_holes_out + k)
                            hole_line_loop_out.append(self.lineId - 1)

                        self.write_line(self.pointId - 2, self.pointId -1 - n_tot)
                        hole_line_loop_out.append(self.lineId - 1)


                    self.write_line_loop(hole_line_loop_out)
                    holes_out.append(self.lineLoopId-1)
                    self.line_loops_out.append(self.lineLoopId - 1)
                    self.write_line_loop(hole_line_loop_in)
                    holes_in.append(self.lineLoopId-1)
                    self.line_loops_in.append(self.lineLoopId - 1)

        self.write_line(self.pointId - 1 - temp, first_point_id)
        line_loop_in_new.append(self.lineId - 1)
        self.write_line_loop(line_loop_in_new)
        self.line_loop_inner_rect.append(self.lineLoopId - 1)
        return holes_in, holes_out

        

    def write_boundary_out(self, topology):
        """
        Here we write the outer rectangle, therfore we write two line loops.
        Holes on the boundary are NOT defined "through" it's own line loop.
        """
        line_loop_in_new = []
        counter = 0
        temp = 0
        t = 0
        for i, edge in enumerate(topology.edges_out()):
            # We write the boundary in following order:
            # left -> up -> right -> down

            for j, hole in enumerate(edge):
                """
                hole is a list.
                if len(list) == 1 than it is just a vertice of the inner rectangle,
                otherwise it is a hole.
                """
                if len(hole) == 1:
                    self.write_point(hole[0][0], hole[0][1], "lc_in")
                    counter += 1

                    if j > 0 or i > 0:
                        self.write_line(self.pointId - temp - 2, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                        temp = 0
                else:
                    s0, s1 = hole[0], hole[-1] # start-  and endpoint of the hole (intersections)
                    hole_in = hole[1:-1]
                    hole_line_loop_in = []

                    # write first intersection point
                    self.write_point(s0[0], s0[1], "lc_in")
                    counter += 1

                    # write connecting line with previous hole
                    if i == 1 and j == 0 and len(topology.edges_out()[0]) != 0:
                        self.write_line(self.pointId - 2 - t, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    elif i > 1 and j == 0:
                        # TODO check bug !

                        self.write_line(self.pointId - 2 - t, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    elif j > 0:
                        self.write_line(self.pointId - 2 - t, self.pointId - 1)
                        line_loop_in_new.append(self.lineId - 1)

                    #write points
                    for h in hole_in:
                        self.write_point(h[0], h[1], "lc_in")
                        counter += 1
                    self.write_point(s1[0], s1[1], "lc_in")
                    counter += 1
                    # write lines
                    n_holes_in = len(hole_in)
                    n_tot = n_holes_in + 2

                    # write lines in
                    for k in range(n_holes_in + 1):
                        self.write_line(self.pointId - n_tot + k, self.pointId - n_tot + k + 1)
                        hole_line_loop_in.append(self.lineId - 1)

                    if s0[i % 2] == s1[i % 2]:
                        # case: edge
                        self.write_line(self.pointId - 1, self.pointId - n_tot)
                        line_loop_in_new.append(- (self.lineId - 1))
                        hole_line_loop_in.append(self.lineId - 1)

                        t = 0
                    else:
                        # case: corner
                        corner_points = [
                        (topology.x0_out, topology.y0_out),
                        (topology.x0_out, topology.y1_out),
                        (topology.x1_out, topology.y1_out),
                        (topology.x1_out, topology.y0_out),
                        ]
                        p = corner_points[i]
                        self.write_point(p[0], p[1], "lc_in")
                        counter += 1
                        temp += 1
                        t = 1
                        # write lines between hole_in and hole_out
                        self.write_line(self.pointId - 1 -1, self.pointId - 1)

                        hole_line_loop_in.append(self.lineId - 1)

                        self.write_line(self.pointId - 1, self.pointId - 1 - n_tot)
                        line_loop_in_new.append(- (self.lineId - 1))
                        line_loop_in_new.append(- (self.lineId - 2))
                        hole_line_loop_in.append(self.lineId - 1)


                    self.write_line_loop(hole_line_loop_in)
                    self.line_loops_out.append(self.lineLoopId - 1)

        
        self.write_line(self.pointId - 1 - t, self.pointId - counter)
        line_loop_in_new.append(self.lineId - 1)
        self.comment("=====================================")
        self.comment("outer line loop:")
        self.comment("=====================================")
        self.write_line_loop(line_loop_in_new)
        self.line_loop_outer_rect = self.lineLoopId - 1
