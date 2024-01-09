import os
import random as rd
from pathlib import Path

# from dolfin import *
import matplotlib.pyplot as plt

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
    topo.mesh(
        file_name=DATA_DIR.joinpath("domain_without_holes_filled"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=True,
    )

    topo.mesh(
        file_name=DATA_DIR.joinpath("domain_without_holes_empty"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=False,
    )


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
    topo.mesh(
        file_name=DATA_DIR.joinpath("2d_1"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=False,
    )


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

    topo.mesh(
        file_name=DATA_DIR.joinpath("2d_3"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=True,
    )


def test_2d_4():
    """Add specific circles."""
    rect_out = [(0, 0), (2, 2)]

    rect_in_1 = [(0.25, 0.25), (0.75, 0.75)]
    rect_in_2 = [(1.25, 0.25), (1.75, 0.75)]
    rect_in_3 = [(1.25, 1.25), (1.75, 1.75)]
    rect_in_4 = [(0.25, 1.25), (0.75, 1.75)]

    topo = Topology(
        rect_out, [rect_in_1, rect_in_2, rect_in_3, rect_in_4], periodic_boundary=False
    )

    radius = 0.05
    midpoints = [
        (1, 0.01),
        (2.01, 0.01),
        (0.251, 0.251),
        (0.251, 0.51),
        (1.251, 1.251),
        (0.4, 0.4),
        (1.751, 1.251),
        (1.251, 1.751),
    ]
    for midpoint in midpoints:
        topo.add_hole(Circle(midpoint, radius), refs=10)

    topo.mesh(
        file_name=DATA_DIR.joinpath("2d_4"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=True,
    )


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

    topo.mesh(
        file_name=DATA_DIR.joinpath("2d_5"),
        write_geo=True,
        lc=1,
        lc_subrects=1,
        filled=True,
    )


def test_2d_6():
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
    N = 10
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

    topo.write_geo(file_name=DATA_DIR.joinpath("2d_6_0"), lc_in=1, lc_out=1, filled=True)
    topo.write_geo(file_name=DATA_DIR.joinpath("2d_6_1"), lc_in=1, lc_out=1, filled=False)
