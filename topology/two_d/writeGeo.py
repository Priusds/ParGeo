# TODO: Check if deprecated.
class WriteGeo(object):
    """
    In this class are basic methods for creating and writing .geo files.
    """
    def __init__(self, file_name, dir = None):
        if dir == None:
            self.file = open(file_name + ".geo", 'w')
        else:
            self.file = open(dir + file_name + ".geo", 'w')

        self.pointId = 1

        self.lineId = 1

        self.lineLoopId = 1

        self.planeSurfaceId = 1

    def write(self,txt):
        self.file.write(txt + ";\n")

    def close(self):
        self.file.close()

    def comment(self, txt):
        self.write("//  " + txt)

    def write_point(self, x, y, lc):
        self.file.write("Point(" + str(self.pointId) + ") = {" + str(x) + "," + str(y) + ",0," + str(
            lc) + "};\n")
        self.pointId += 1

    def write_line(self, p1, p2):
        self.file.write("Line(" + str(self.lineId) + ") = {" + str(p1) + "," + str(p2) + "};\n")
        self.lineId += 1

    def write_line_loop(self, line_loop):
        self.file.write("Line Loop(" + str(self.lineLoopId) + ") = {" + ','.join(
            [str(i) for i in line_loop]) + "};\n")
        self.lineLoopId += 1

    def write_plane_surface(self, loops):
        self.file.write("Plane Surface(" + str(self.planeSurfaceId) + ") = {" + ','.join(
            [str(i) for i in loops]) + "};\n")
        self.planeSurfaceId += 1

    def write_polygon(self, polygon, lc):
        n = len(polygon)
        # write all points
        for vertex in polygon:
            x, y = vertex
            self.write_point(x, y, lc)
        # write all lines
        for i in range(n - 1):
            self.write_line(self.pointId - n + i, self.pointId - n + i + 1)
        self.write_line(self.pointId - 1, self.pointId - n)
        # write line loop
        self.write_line_loop(range(self.lineId - n, self.lineId))

    def write_physical_surface(self, tag, surface_Id):
        if isinstance(surface_Id, str):
            self.file.write("Physical Surface(" + tag + ") = {" + surface_Id + "};\n")
        else:
            if len(surface_Id) != 0:
                self.file.write("Physical Surface(" + tag + ") = {" + ','.join(surface_Id) + "};\n")
            else:
                pass
