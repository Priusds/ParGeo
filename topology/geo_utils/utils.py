def write_point(data):
    lc = 1 if data["lc"] is None else data["lc"]
    string = "Point({0}) = {{{1},{2},{3},{4}}};\n".format(data["id"],
                                            data["coords"][0],
                                            data["coords"][1],
                                            data["coords"][2],
                                            lc)
    return string

def write_line(data):
    return "Line({0}) = {{{1},{2}}};\n".format(data["id"],data["start"],data["end"])

def write_lineLoop(data):
    loop = ""
    for i in data["lines"]:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Line Loop({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_surface(data):
    loop = ""
    for i in data['lineLoops']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Plane Surface({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_surfaceLoop(data):
    loop = ""
    for i in data['surfaces']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Surface Loop({0}) = {{{1}}};\n".format(data["id"], loop)
    return string

def write_volume(data):
    loop = ""
    for i in data['surfaceLoops']:
        if len(loop)>0:
            loop += ","
        loop += str(i)
    string = "Volume({0}) = {{{1}}};\n".format(data["id"], loop)
    return string


functions = {
    "points":write_point,
    "lines":write_line,
    "lineLoops":write_lineLoop,
    "surfaces":write_surface,
    "surfaceLoops":write_surfaceLoop,
    "volumes":write_volume
}

def write_geo(geo, file_name):
    file = open(file_name + ".geo", 'w')
    for key in ["points", "lines","lineLoops", "surfaces", "surfaceLoops","volumes"]: #,"surfaces"]:
        if key in geo:
            for value in geo[key]:
                file.write(functions[key](value))
    file.close()