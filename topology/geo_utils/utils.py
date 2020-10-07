def write_point(Id, data):
    lc = 1 if data["lc"] is None else data["lc"]
    string = "Point({0}) = {{{1},{2},{3},{4}}};\n".format(Id,
                                            data["coords"][0],
                                            data["coords"][1],
                                            data["coords"][2],
                                            lc)
    return string

def write_line(Id, data):
    return "Line({0}) = {{{1},{2}}};\n".format(Id,data[0],data[1])

def write_lineLoop(Id, data):
    loop = ""
    for i in data:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Line Loop({0}) = {{{1}}};\n".format(Id, loop)
    return string

def write_surface(Id, data):
    loop = ""
    for i in data:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Plane Surface({0}) = {{{1}}};\n".format(Id, loop)
    return string

def write_physicalSurface(Id, data):
    loop = ""
    for i in data:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Physical Surface({0}) = {{{1}}};\n".format(Id, loop)
    return string

def write_surfaceLoop(Id, data):
    loop = ""
    for i in data:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Surface Loop({0}) = {{{1}}};\n".format(Id, loop)
    return string

def write_volume(Id, data):
    loop = ""
    for i in data:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Volume({0}) = {{{1}}};\n".format(Id, loop)
    return string


functions = {
    "points":write_point,
    "lines":write_line,
    "lineLoops":write_lineLoop,
    "surfaces":write_surface,
    "physicalSurfaces":write_physicalSurface,
    "surfaceLoops":write_surfaceLoop,
    "volumes":write_volume
}

def write_geo(geo, file_name):
    file = open(file_name + ".geo", 'w')
    for key in ["points", "lines","lineLoops", "surfaces", "surfaceLoops","volumes", "physicalSurfaces"]: 
        if key in geo:
            for Id, value in geo[key].items():
                file.write(functions[key](Id, value))
    file.close()

def write_geometry(geometry, file_name):
    file = open(file_name + ".geo", 'w')
    
##########################

def write_point2(data):
    lc = 1 if data["lc"] is None else data["lc"]
    string = "Point({0}) = {{{1},{2},{3},{4}}};\n".format(data["id"],
                                            data["coords"][0],
                                            data["coords"][1],
                                            data["coords"][2],
                                            lc)
    return string

def write_line2(data):
    return "Line({0}) = {{{1},{2}}};\n".format(data["id"],data["start"],data["end"])

def write_lineLoop2(data):
    loop = ""
    for i in data["lines"]:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Line Loop({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_surface2(data):
    loop = ""
    for i in data['lineLoops']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Plane Surface({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_surfaceLoop2(data):
    loop = ""
    for i in data['surfaces']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Surface Loop({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_volume2(data):
    loop = ""
    for i in data['surfaceLoops']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Volume({0}) = {{{1}}};\n".format(data["id"], loop)
    return string



def write_geo2(geo, file_name):
    file = open(file_name + ".geo", 'w')
    for key in ["points", "lines","lineLoops", "surfaces", "surfaceLoops","volumes"]: #,"surfaces"]:
        if key in geo:
            for value in geo[key]:
                file.write(functions[key](value))
    file.close()

def write_geometry2(geometry, file_name):
    file = open(file_name + ".geo", 'w')
    