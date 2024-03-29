"""Move all dolfine related code from Topology to here."""
import os
import dolfin

def geo_to_xml(file_name, width, height):
    os.system("gmsh -2 " + file_name + ".geo")  # -v 0
    os.system("dolfin-convert " + file_name + ".msh " + file_name + ".xml")
    mesh = dolfin.Mesh(file_name + ".xml")
    sort_xml(file_name, mesh, width, height)

def sort_xml(file_name, mesh, a, b):

    ct = {} # coordinate transformation for the boundary

    left_y  = sorted([p for p in mesh.coordinates() if p[0] == 0], key = lambda x: x[1])
    right_y = sorted([p for p in mesh.coordinates() if p[0] == b], key = lambda x: x[1])
    up_x    = sorted([p for p in mesh.coordinates() if p[1] == a], key = lambda x: x[0])
    down_x  = sorted([p for p in mesh.coordinates() if p[1] == 0], key = lambda x: x[0])

    assert len(left_y) == len(right_y) and len(up_x) == len(down_x)

    l_y = len(left_y)
    l_x = len(down_x)

    for p_l, p_r in zip(left_y, right_y):
        # construct old point pair as key ...
        keys = [(p_l[0], p_l[1]), (p_r[0], p_r[1])]
        new_y = 0.5*(p_l[1]+p_r[1])
        for key in keys:
            # and map it to the pertubated new point
            ct[key] = (key[0], new_y)


    for p_b, p_t in zip(down_x, up_x):
        # construct old point pair as key ...
        keys = [(p_b[0], p_b[1]), (p_t[0], p_t[1])]
        new_x = 0.5*(p_b[0]+p_t[0])
        for key in keys:
            # and map it to the pertubated new point
            ct[key] = (new_x, key[1])

    for coord in mesh.coordinates():

        key = (coord[0],coord[1])
        x = coord[0]
        y = coord[1]
        if x == 0 or x == a or y == 0 or y == b:

            assert(key in ct)
            new_coord = ct[(x,y)]
            coord[0] = new_coord[0]
            coord[1] = new_coord[1]


    # now save the mesh to xml ...
    File(file_name + ".xml") << mesh
