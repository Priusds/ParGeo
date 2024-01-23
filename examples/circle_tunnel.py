import shapely
from bubbles.gmsh_api import write_geo, topology_to_gmsh_entities
from bubbles.two_d.geometry import Circle
from bubbles.two_d.topology import Topology

def generate_topo():
    """Generate a topology with a sequence of centered inclusions."""
    C = shapely.Polygon(Circle(midpoint=(0., 0), radius=1.).discretize_hole(refs=50))

    topo = Topology(C)

    radii   = [0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
    level = [i+1 for i in range(len(radii))]

    for lvl, rad in zip(level, radii): 
        C = shapely.Polygon(Circle(midpoint=(0, 0), radius=rad).discretize_hole(refs=50))
        topo.add(polygon=C, level= lvl)

    return topo

if __name__ == "__main__":
    topo = generate_topo()
    topo.plot()

    gmsh_entities = topology_to_gmsh_entities(topo)
    write_geo(
        gmsh_entities=gmsh_entities, file_name="circle_tunnel", correct_curve_loops=True
    )
