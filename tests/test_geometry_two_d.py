import os
import random as rd
from pathlib import Path

# from dolfin import *
import matplotlib.pyplot as plt

from bubbles.geo_utils.utils import write_geo
from bubbles.two_d.hole import Circle, Stellar
from bubbles.two_d.topology import Topology

# Safe .GEO files at DATA_DIR
DATA_DIR = Path(__file__).parent.joinpath("data/geometry_two_d")
DATA_DIR.mkdir(exist_ok=True, parents=True)

"""TODO: Make dolfin dependent test optional
def dolfin_mesh(file_name):
    # creates .xml mesh from .geo mesh and returns dolfin mesh
    os.system("gmsh -2 " + file_name + ".geo")  # -v 0
    os.system("dolfin-convert " + file_name + ".msh " + file_name + ".xml")
    return Mesh(file_name + ".xml")
"""


def test_domain_without_holes_fill():
    """Create a domain without holes and filled sub-rects."""
    rect_out = [(0, 0), (2, 2)]
    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=True
    )
    geo = topo.get_geometry(lc=1.0, lc_subrects=0.1, filled=True)
    write_geo(geo, DATA_DIR / "domain_without_holes_filled")


def test_domain_without_holes_empty():
    """Create a domain without holes and empty sub-rects."""
    rect_out = [(0, 0), (2, 2)]
    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=True
    )
    geo = topo.get_geometry(lc=1.0, lc_subrects=0.1, filled=False)
    write_geo(geo, DATA_DIR / "domain_without_holes_empty")


def test_2d_1():
    """Domain with four sub-rects and many stellar holes."""
    # Domain with many little stellar holes
    # Define the rectangles
    width = 2
    height = 2
    rect_out = [(0, 0), (2, 2)]
    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    # Set periodic_boundary, True or False for outer rectangle

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=True
    )

    # Make some holes
    N = 100
    counter = 0
    i = 0
    while counter < N and i < 1000:
        i += 1
        print(i)
        midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
        # midpoint_s = (.5,1) # TODO: check error!
        # this way you define a Stellar hole
        stellar_hole = Stellar(midpoint_s, 0.02)

        if topo.add_hole(stellar_hole, refs=30):
            counter += 1

    print("Needed ", i, "runs for ", N, "holes")

    # Make geometry, set filled True or False if the inclusion is a hole or not
    geo = topo.get_geometry(lc=1, lc_subrects=1, filled=False)
    write_geo(geo, DATA_DIR.joinpath("2d_1"))


"""TODO: Make dolfin dependent test optional
def test_2d_2():
    file_name = "test"
    mesh = dolfin_mesh("test")
    subdomains_hole = MeshFunction("size_t", mesh, "test" + "_physical_region.xml")
    plot(subdomains_hole)
    #plot(mesh)
    plt.show()
"""


def test_2d_3():
    """Domain with one sub-rect and two specific holes."""
    rect_out = [(0, 0), (2, 2)]
    rect_in = [(0.5, 0.5), (1.5, 1.5)]

    topo = Topology(rect_out, [rect_in])

    midpoint_s = (0.5, 1)
    stellar_hole = Circle(midpoint_s, 0.05)
    topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (0.51, 1.51)
    stellar_hole = Circle(midpoint_s, 0.05)
    topo.add_hole(stellar_hole, refs=10)

    topo.write_geo(DATA_DIR.joinpath("2d_3"))


def test_2d_4():
    rect_out = [(0, 0), (2, 2)]

    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=False
    )

    midpoint_s = (1, 0.01)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (2.01, 0.01)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=20)

    midpoint_s = (0.251, 0.251)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (0.251, 0.51)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (1.251, 1.251)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (0.4, 0.4)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (1.751, 1.251)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    midpoint_s = (1.251, 1.751)
    stellar_hole = Circle(midpoint_s, 0.05)
    added = topo.add_hole(stellar_hole, refs=10)

    geo = topo.get_geometry(lc=1, lc_subrects=1, filled=True)
    write_geo(geo, DATA_DIR.joinpath("2d_4"))


def test_2d_5():
    # Domain with many little stellar holes
    # Define the rectangles
    width = 2
    height = 2
    rect_out = [(0, 0), (2, 2)]
    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    # Set periodic_boundary, True or False for outer rectangle

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=True
    )

    # Make some holes
    N = 50
    counter = 0
    i = 0
    while counter < N and i < 1000:
        i += 1
        print(i)
        midpoint_s = rd.uniform(0, width), rd.uniform(0, height)
        radius = rd.uniform(0, 0.5)
        stellar_hole = Circle(midpoint_s, radius)

        if topo.add_hole(stellar_hole, refs=30):
            counter += 1

    print("Needed ", i, "runs for ", N, "holes")

    # Make geometry, set filled True or False if the inclusion is a hole or not
    geo = topo.get_geometry(lc=1, lc_subrects=1, filled=True)
    write_geo(geo, DATA_DIR.joinpath("2d_5"))
