import math
import shapely
from itertools import groupby

from bubbles.gmsh_api import Point, Line, CurveLoop, PlaneSurface, PhysicalGroup, SurfaceLoop, Volume, GmshEntities, mesh, write_geo
from bubbles.two_d.hole import Circle, Ellipse, Rectangular
from bubbles.two_d.topology_v2 import Bubble, Topology

def generate_topo():
    R = Rectangular(midpoint=(0.0, 0), width=2, height=2).discretize_hole(refs=4)
    C = Circle(midpoint=(1.0, 0), radius=0.5).discretize_hole(refs=50)

    C_cutout = shapely.Polygon(Circle(midpoint=(0. ,0.75), radius = 0.2).discretize_hole(refs=50))

    R_shapely = shapely.Polygon(R)
    C_shapely = shapely.Polygon(C)
    domain = R_shapely#.union(C_shapely)#.difference(C_cutout)
    topo = Topology(domain)

    if True:
        # P2
        c1 = shapely.Polygon(
            Circle(midpoint=(-0.25, 0.), radius=0.25).discretize_hole(refs=50)
        )
        topo.add_bubble(c1, level=10, is_hole=True)
    
    if True:
        c2 = shapely.Polygon(Circle(midpoint=(-0.1 ,0.4), radius = 0.3).discretize_hole(refs=50))
        topo.add_bubble(c2, level=11, is_hole=False)

    return topo

    # P3
    C3 = Circle(midpoint=(0.75, 0), radius=0.25).discretize_hole(refs=50)
    p3 = shapely.Polygon(C3)
    topo.add_bubble(p3, level=1, is_hole=False)

    # P4
    C4 = Circle(midpoint=(1.0, 0), radius=0.25).discretize_hole(refs=50)
    p4 = shapely.Polygon(C4)
    topo.add_bubble(p4, level=1, is_hole=False)
    
    # this guy is intersection with the domain
    E1 = Ellipse(
        midpoint=(-0.5, 0.75), axis=(1, 0.1), angle=0.25 * math.pi
    ).discretize_hole(refs=50)
    
    return topo
    # this guy is just interior
    E2 = Ellipse(midpoint=(-0.55, -0.6), axis=(0.25, 0.1), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E1), level=1, is_hole=False)
    topo.add_bubble(shapely.Polygon(E2), level=1, is_hole=False)

    


    C5 = Circle(midpoint=(0.5, 0.6), radius=0.2).discretize_hole(refs=50)
    C6 = Circle(midpoint=(0.6, 0.6), radius=0.2).discretize_hole(refs=50)
    topo.add_bubble(shapely.Polygon(C5), level=3, is_hole=False)
    topo.add_bubble(shapely.Polygon(C6), level=2, is_hole=False)

    C = Circle(midpoint=(0.24, -0.8), radius=0.22).discretize_hole(refs=50)
    topo.add_bubble(shapely.Polygon(C), level=1, is_hole=False)


 


    E = Ellipse(midpoint=(0.25, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)

    E = Ellipse(midpoint=(0.4, -0.4), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)

    E = Ellipse(midpoint=(0.5, -0.6), axis=(0.05, 0.2), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)



    # this guy encloses the loop !
    E = Ellipse(midpoint=(0.4, -0.7), axis=(0.2, 0.05), angle=0.0).discretize_hole(
        refs=50
    )
    topo.add_bubble(shapely.Polygon(E), level=2, is_hole=False)
    return topo


def bubble_to_gmsh_entities(bubble: Bubble, point_tag, line_tag, curve_loop_tag, plane_surface_tag):
    assert not bubble.is_hole
    points_total = []
    lines_total = []
    curve_loops = []


    
    physical_group_tag = bubble.level + 1 # +1 because physical group tags start at 1
    # ==================================================
    # Write exterior loop
    # ==================================================
    assert isinstance(bubble.polygon, shapely.Polygon)
    # Add points
    points_local = []
    for x,y in  bubble.polygon.exterior.coords[:-1]:
        point = Point(x=x, y=y, z=0.0, lc=0.1, tag=point_tag)
        points_local.append(point)
        point_tag += 1

    # Add lines
    lines_local = []
    for i in range(len(points_local)-1):
        line = Line(start_tag=points_local[i].tag, end_tag=points_local[i+1].tag, tag=line_tag)
        lines_local.append(line)
        line_tag += 1
    line = Line(start_tag=points_local[-1].tag, end_tag=points_local[0].tag, tag=line_tag)
    lines_local.append(line)
    line_tag += 1

    points_total = points_total + points_local
    lines_total = lines_total + lines_local
    # Add curve loop
    curve_loop = CurveLoop(line_tags=[line.tag for line in lines_local], tag=curve_loop_tag)
    curve_loops.append(curve_loop)
    curve_loop_tag += 1

    # ==================================================
    # Write interior loops
    # ==================================================
    for loop in bubble.polygon.interiors:
        points_local = []
        lines_local = []
        for x,y in  loop.coords[:-1]:
            point = Point(x=x, y=y, z=0.0, lc=0.1, tag=point_tag)
            points_local.append(point)
            point_tag += 1

        # Add lines
        for i in range(len(points_local)-1):
            line = Line(start_tag=points_local[i].tag, end_tag=points_local[i+1].tag, tag=line_tag)
            lines_local.append(line)
            line_tag += 1
        line = Line(start_tag=points_local[-1].tag, end_tag=points_local[0].tag, tag=line_tag)
        lines_local.append(line)
        line_tag += 1
        
        # Add curve loop
        curve_loop = CurveLoop(line_tags=[line.tag for line in lines_local], tag=curve_loop_tag)
        curve_loops.append(curve_loop)
        curve_loop_tag += 1

        points_total = points_total + points_local
        lines_total = lines_total + lines_local
    
    # ==================================================
    # Write surface
    # ==================================================
    plane_surfaces = [PlaneSurface(curve_loop_tags=[loop.tag for loop in curve_loops], tag=plane_surface_tag)]
    physical_groups = [PhysicalGroup(dim=2, entity_tags=[plane_surface_tag], tag=physical_group_tag)]
    plane_surface_tag += 1

    return GmshEntities(
        points=points_total,
        lines=lines_total,
        curve_loops=curve_loops,
        plane_surfaces=plane_surfaces,
        physical_groups=physical_groups,
    ), point_tag, line_tag,curve_loop_tag, plane_surface_tag

def topology_to_gmsh_entities(topo: Topology):
    all_entities_list = []

    bubbles = list(topo.bubbles)

    point_tag = 1
    line_tag = 1
    curve_loop_tag = 1
    plane_surface_tag = 1

    for bubble in bubbles:
        if not bubble.is_hole:
            gmsh_entities, point_tag, line_tag,curve_loop_tag, plane_surface_tag = bubble_to_gmsh_entities(bubble, point_tag, line_tag, curve_loop_tag, plane_surface_tag)
            all_entities_list.append(gmsh_entities)
    
    # Group all physical groups by tag
    physical_groups_merged = []
    physical_groups = [e.physical_groups[0] for e in all_entities_list]
    physical_groups.sort(key=lambda x: x.tag)


    for tag, physical_group in groupby(physical_groups, lambda x: x.tag):
        entity_tags = [s.entity_tags[0] for s in physical_group]
        physical_groups_merged.append(PhysicalGroup(dim=2, entity_tags=entity_tags, tag=tag))

    return GmshEntities(
        points=[i for ent in all_entities_list for i in ent.points],
        lines=[i for ent in all_entities_list for i in ent.lines],
        curve_loops=[i for ent in all_entities_list for i in ent.curve_loops],
        plane_surfaces=[i for ent in all_entities_list for i in ent.plane_surfaces],
        physical_groups=physical_groups_merged,
    )


if __name__ == "__main__":
    for _ in range(1):
        generate_topo()
    #topo = generate_topo()
    #topo.plot_setup()

    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(gmsh_entities=gmsh_entities, file_name="mesh", correct_curve_loops=True)
    quit()
    mesh(
        gmsh_entities=gmsh_entities_all,
        dim=2,
        file_name="mesh",
        save_all=False,
        write_geo=True,
        correct_curve_loops=True,
    )