
import random

import shapely
from bubbles.constraint import DistanceConstraint
from bubbles.geometry import Circle, Rectangle, Stellar, Ellipse
from bubbles.gmsh_utils import topology_to_gmsh_entities, write_geo
from bubbles.topology import Topology

import math


def test():
        __alpha = math.pi / 6
        domain = Ellipse(midpoint=(0., 0.),axis = (0.5,1), angle = 0.3).discretize(refs=50)
        topology = Topology(domain)

        #polygon =  Stellar(midpoint=(0,0), radius=0.025).discretize(refs=50)
        polygon =  Ellipse(midpoint=(0,0), axis=(0.01,0.04), angle = 0).discretize(refs=50)
        #Rectangle(midpoint=(0., 0.), width=0.05, height=0.05).to_polygon()
        
        Lx = 0.05
        Ly = 0.05
      

       # center of bounding box is rotational origin
        mx, my = -0.5, 0.3

        rotate_x = lambda x,y : math.cos(__alpha) *x - math.sin(__alpha) *y
        rotate_y = lambda x,y : math.sin(__alpha) *x + math.cos(__alpha) *y

        domain_rotated = shapely.Polygon([(rotate_x(x-mx,y-my)+mx,rotate_y(x-mx,y-my)+my) for x,y in zip(*domain.exterior.xy)])
        polygon_rotated = shapely.Polygon([(rotate_x(x-mx,y-my)+mx,rotate_y(x-mx,y-my)+my) for x,y in zip(*polygon.exterior.xy)])

        #x,y = domain_rotated.exterior.xy

        enclosing_circle = shapely.minimum_rotated_rectangle(domain)
        minx, miny, maxx, maxy  = enclosing_circle.bounds

        #minx, miny, maxx, maxy = domain_rotated.bounds
        poly_minx, poly_miny, poly_maxx, poly_maxy = polygon_rotated.bounds

        bounding_box_rotated_poly = shapely.box(poly_minx, poly_miny, poly_maxx, poly_maxy)
        bounding_box_rotated_poly = shapely.Polygon([ (rotate_x(x-mx,y-my)+mx,rotate_y(x-mx,y-my)+my) for x,y in zip(*bounding_box_rotated_poly.exterior.xy) ])


        # poly_minx + kx_max * Lx <= maxx <=>   (maxx - poly_minx)/Lx >= kx_max
        kx_max = math.floor((maxx - poly_minx) / Lx)
        # poly_maxx + kx_min * Lx >= minx  <=>  (minx - poly_maxx)/Lx <= kx_min
        kx_min = math.ceil((minx - poly_maxx) / Lx)
        # analog for y axis
        ky_max = math.floor((maxy - poly_miny) / Ly)
        ky_min = math.ceil((miny - poly_maxy) / Ly)

        # collect possible polygons for periodic rule
        # periodic_polygons = [
        #     shapely.affinity.translate(polygon, xoff=rotate_x(kx * Lx-mx, ky * Ly -my)+mx, yoff= rotate_y(kx * Lx-mx,ky * Ly -my)+my)
        #     for kx in range(kx_min, kx_max + 1)
        #     for ky in range(ky_min, ky_max + 1)
        # ]

        bounding_poly = shapely.box(poly_minx, poly_miny, poly_maxx, poly_maxy)
        periodic_polygons =[]
        for kx in range(kx_min, kx_max + 1):
            for ky in range(ky_min, ky_max + 1):
                   box_rotate = shapely.affinity.translate(bounding_poly, xoff= rotate_x(kx * Lx-mx, ky * Ly -my) + mx,
                                                                          yoff= rotate_y(kx * Lx-mx, ky * Ly -my) + my)
                   if box_rotate.intersects(domain):
                        periodic_polygons.append(shapely.affinity.translate(polygon, xoff= rotate_x(kx * Lx-mx, ky * Ly -my)+mx,
                                                                                     yoff= rotate_y(kx * Lx-mx, ky * Ly -my)+my))


        # # and possible union them
        # periodic_polygons = shapely.union_all(
        #     periodic_polygons, grid_size=topology.grid_size
        # )

        import matplotlib.pyplot as plt
        x,y = domain.exterior.xy
        plt.plot(x,y, '-k')

        for poly in periodic_polygons:
            x,y = poly.exterior.xy
            plt.plot(x,y, '-b')

        x,y = bounding_box_rotated_poly.exterior.xy
        plt.plot(x,y, '-r')

       

        x,y = enclosing_circle.exterior.xy
        plt.plot(x,y, '-r')
              
        plt.show()



if __name__ == "__main__":
    test()
    #topo = generate_topo_football()  # generate_topo_simple() #generate_topo()
    #topo.plot()
    # gmsh_entities = topology_to_gmsh_entities(topo)
    # write_geo(
    #     gmsh_entities=gmsh_entities,
    #     file_name="stellar_samples_periodic_constrains",
    #     correct_curve_loops=True,
    # )
