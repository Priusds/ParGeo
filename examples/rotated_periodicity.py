import math

import shapely

from pargeo.geometry import Ellipse, Stellar
from pargeo.topology import Topology


def test():
    __alpha = math.pi / 6
    domain = Ellipse(midpoint=(0.0, 0.0), axis=(0.5, 1), angle=0.3).discretize(refs=50)


    # polygon =  Stellar(midpoint=(0,0), radius=0.025).discretize(refs=50)
    polygon = Ellipse(midpoint=(0, 0), axis=(0.01, 0.04), angle=0).discretize(refs=50)
    # Rectangle(midpoint=(0., 0.), width=0.05, height=0.05).to_polygon()

    Lx = 0.05
    Ly = 0.05

    # center of bounding box is rotational origin
    mx, my = -0.5, 0.3

    rotate_x = lambda x, y: math.cos(__alpha) * x - math.sin(__alpha) * y
    rotate_y = lambda x, y: math.sin(__alpha) * x + math.cos(__alpha) * y

    domain_rotated = shapely.Polygon(
        [
            (rotate_x(x - mx, y - my) + mx, rotate_y(x - mx, y - my) + my)
            for x, y in zip(*domain.exterior.xy)
        ]
    )
    polygon_rotated = shapely.Polygon(
        [
            (rotate_x(x - mx, y - my) + mx, rotate_y(x - mx, y - my) + my)
            for x, y in zip(*polygon.exterior.xy)
        ]
    )

    # x,y = domain_rotated.exterior.xy

    enclosing_circle = shapely.minimum_rotated_rectangle(domain)
    minx, miny, maxx, maxy = enclosing_circle.bounds

    # minx, miny, maxx, maxy = domain_rotated.bounds
    poly_minx, poly_miny, poly_maxx, poly_maxy = polygon_rotated.bounds

    bounding_box_rotated_poly = shapely.box(poly_minx, poly_miny, poly_maxx, poly_maxy)
    bounding_box_rotated_poly = shapely.Polygon(
        [
            (rotate_x(x - mx, y - my) + mx, rotate_y(x - mx, y - my) + my)
            for x, y in zip(*bounding_box_rotated_poly.exterior.xy)
        ]
    )

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
    periodic_polygons = []
    for kx in range(kx_min, kx_max + 1):
        for ky in range(ky_min, ky_max + 1):
            box_rotate = shapely.affinity.translate(
                bounding_poly,
                xoff=rotate_x(kx * Lx - mx, ky * Ly - my) + mx,
                yoff=rotate_y(kx * Lx - mx, ky * Ly - my) + my,
            )
            if box_rotate.intersects(domain):
                periodic_polygons.append(
                    shapely.affinity.translate(
                        polygon,
                        xoff=rotate_x(kx * Lx - mx, ky * Ly - my) + mx,
                        yoff=rotate_y(kx * Lx - mx, ky * Ly - my) + my,
                    )
                )

    # # and possible union them
    # periodic_polygons = shapely.union_all(
    #     periodic_polygons, grid_size=topology.grid_size
    # )

    import matplotlib.pyplot as plt

    x, y = domain.exterior.xy
    plt.plot(x, y, "-k")

    for poly in periodic_polygons:
        x, y = poly.exterior.xy
        plt.plot(x, y, "-b")

    x, y = bounding_box_rotated_poly.exterior.xy
    plt.plot(x, y, "-r")

    x, y = enclosing_circle.exterior.xy
    plt.plot(x, y, "-r")

    plt.show()


if __name__ == "__main__":
    test()
    # topo = generate_topo_football()  # generate_topo_simple() #generate_topo()
    # topo.plot()
    # gmsh_entities = topology_to_gmsh_entities(topo)
    # write_geo(
    #     gmsh_entities=gmsh_entities,
    #     file_name="stellar_samples_periodic_constrains",
    #     correct_curve_loops=True,
    # )

import math

import shapely

from pargeo.geometry import Ellipse, Stellar
from pargeo.topology import Topology


def test():
    v1, v2 = 1.0, 0.0
    w1, w2 = 2.0, 1.0

    v1, v2 = v1 / math.sqrt(v1**2 + v2**2), v2 / math.sqrt(v1**2 + v2**2)
    w1, w2 = w1 / math.sqrt(w1**2 + w2**2), w2 / math.sqrt(w1**2 + w2**2)

    # __alpha = math.pi / 6
    Lv = 0.1
    Lw = 0.1

    # domain = Ellipse(midpoint=(0., 0.),axis = (0.5,1), angle = 0.3).discretize(refs=100)

    domain = Stellar(midpoint=(0.0, 0.0), radius=1).discretize(refs=100)
    topo = Topology(domain)

    # polygon =  Stellar(midpoint=(0,0), radius=0.025).discretize(refs=50)
    polygon = Ellipse(midpoint=(0, 0), axis=(0.01, 0.01), angle=-0.5).discretize(
        refs=50
    )
    # Rectangle(midpoint=(0., 0.), width=0.05, height=0.05).to_polygon()

    d = domain.bounds[2] - domain.bounds[0] + domain.bounds[3] - domain.bounds[1]
    Kv = math.ceil(d / Lv)
    Kw = math.ceil(d / Lw)

    # r_dom = shapely.minimum_bounding_radius(domain)
    # r_poly = shapely.minimum_bounding_radius(polygon)
    # v = math.ceil(2*r_dom / Lv)
    # d = domain.bounds[2] - domain.bounds[0] + domain.bounds[3] - domain.bounds[1]

    v_dirs = [(v1, v2), (-v1, -v2)]
    w_dirs = [(w1, w2), (-w1, -w2)]

    # v_dirs = [(v1,v2)]
    # w_dirs = [(w1,w2)]

    level = 1
    periodic_polygons = shapely.union_all(
        [
            Repeat(v_dir, Lv, Kv, w_dir, Lw, Kw)(polygon, level, topo)
            for v_dir in v_dirs
            for w_dir in w_dirs
        ]
    )

    # periodic_polygons = []

    # K = max(Kv,Kw)

    # offsets = [ (i * Lv * v1 + j * Lw * w1,
    #              i * Lv * v2 + j * Lw * w2)  for i in range(-Kv, Kv+1) for j in range(-Kw, Kw+1)
    #              if (i * Lv * v1 + j * Lw * w1)**2 + (i * Lv * v2 + j * Lw * w2)**2 <= (2*r_dom)**2 + r_poly**2]
    # for x_off,y_off in offsets:
    #     poly_shifted = shapely.affinity.translate(polygon, xoff =x_off, yoff =  y_off)
    #     if poly_shifted.intersects(domain):
    #             periodic_polygons.append(poly_shifted)

    # for kv in range(-Kv, Kv+1):
    #      for kw in range(-Kw, Kw+1):
    # # for kv in range(-K, K+1):
    # #      for kw in range(-K, K+1):
    #         # kv*Lv * (v1,v2)  + kw * Lw * (w1,w2)
    #         x_off = kv * Lv * v1 + kw * Lw * w1
    #         y_off = kv * Lv * v2 + kw * Lw * w2
    #         #if x_off**2 + y_off**2 > (2*r_dom)**2 + r_poly**2: # TODO: is the 2r_dom necessary?
    #         #     continue
    #         #else:
    #         poly_shifted = shapely.affinity.translate(polygon, xoff =x_off, yoff =  y_off)
    #         if poly_shifted.intersects(domain):
    #                 periodic_polygons.append(poly_shifted)

    # periodic_polygons = shapely.union_all([ Repeat( (d1,Lv,Kv), (d2,Lw,Kw))(polygon, level, topo) for d1,d2 in directions])

    # for kv in range(0, )
    # for kv, kw in index_set():
    #     box_rotate = shapely.affinity.translate(enclosing_circle, xoff = kv * Lv * v1 + kw * Lw * w1 ,
    #                                                               yoff =  kv * Lv * v2 + kw * Lw * w2)
    #     if box_rotate.intersects(domain):
    #             periodic_polygons.append(shapely.affinity.translate(polygon, xoff= kv * Lv * v1 + kw * Lw * w1,
    #                                                                          yoff= kv * Lv * v2 + kw * Lw * w2))

    # #minx, miny, maxx, maxy = domain_rotated.bounds
    # poly_minx, poly_miny, poly_maxx, poly_maxy = polygon_rotated.bounds

    # bounding_box_rotated_poly = shapely.box(poly_minx, poly_miny, poly_maxx, poly_maxy)
    # bounding_box_rotated_poly = shapely.Polygon([ (rotate_x(x-mx,y-my)+mx,rotate_y(x-mx,y-my)+my) for x,y in zip(*bounding_box_rotated_poly.exterior.xy) ])

    # poly_minx + kx_max * Lx <= maxx <=>   (maxx - poly_minx)/Lx >= kx_max
    # kx_max = math.floor((maxx - poly_minx) / Lv)
    # # poly_maxx + kx_min * Lx >= minx  <=>  (minx - poly_maxx)/Lx <= kx_min
    # kx_min = math.ceil((minx - poly_maxx) / Lv)
    # # analog for y axis
    # ky_max = math.floor((maxy - poly_miny) / Lw)
    # ky_min = math.ceil((miny - poly_maxy) / Lw)

    # collect possible polygons for periodic rule
    # periodic_polygons = [
    #     shapely.affinity.translate(polygon, xoff=rotate_x(kx * Lx-mx, ky * Ly -my)+mx, yoff= rotate_y(kx * Lx-mx,ky * Ly -my)+my)
    #     for kx in range(kx_min, kx_max + 1)
    #     for ky in range(ky_min, ky_max + 1)
    # ]

    # bounding_poly = shapely.box(poly_minx, poly_miny, poly_maxx, poly_maxy)
    # periodic_polygons =[]
    # for kx in range(kx_min, kx_max + 1):
    #     for ky in range(ky_min, ky_max + 1):
    #            box_rotate = shapely.affinity.translate(bounding_poly, xoff= rotate_x(kx * Lx-mx, ky * Ly -my) + mx,
    #                                                                   yoff= rotate_y(kx * Lx-mx, ky * Ly -my) + my)
    #            if box_rotate.intersects(domain):
    #                 periodic_polygons.append(shapely.affinity.translate(polygon, xoff= rotate_x(kx * Lx-mx, ky * Ly -my)+mx,
    #                                                                              yoff= rotate_y(kx * Lx-mx, ky * Ly -my)+my))

    # # and possible union them
    # periodic_polygons = shapely.union_all(
    #     periodic_polygons, grid_size=topology.grid_size
    # )

    import matplotlib.pyplot as plt

    x, y = domain.exterior.xy
    plt.plot(x, y, "-k")

    for poly in periodic_polygons.geoms:
        x, y = poly.exterior.xy
        plt.plot(x, y, "-b")

    # x,y = bounding_box_rotated_poly.exterior.xy
    # plt.plot(x,y, '-r')

    # x,y = enclosing_circle.exterior.xy
    # plt.plot(x,y, '-r')

    plt.show()


if __name__ == "__main__":
    test()
