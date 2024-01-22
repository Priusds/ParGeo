import shapely


class Distance_Constraint(object):
    def set_level_distance(self, distance, lvl_set_1, lvl_set_2):
        # lvl_1 = int | list[int]
        # lvl_2 = int | list[int]
        # distance = float | list[float]

        # logic: if distance = float apply distance for all (symmetric) combinations of lvl_1, lvl_2
        #        if distance = list, len must match  lvl_1 or lvl_2  list len

        # store a (symmetric) hashmap of distance values that can be overwritten with warning !

        pass

    def set_level_all(self, infos: dict()):
        pass

    def set_geometry_distance(self, level: shapely.Geometry, distance):
        pass

    def __call__(self, polygon):
        pass


# TODO: Just add a collision: poly -> bool argument to the add function, it could for example
# check the number of refs.
def polygon_distance_req(
    poly1: shapely.Polygon | shapely.MultiPolygon,
    poly2: shapely.Polygon | shapely.MultiPolygon,
    distance: float,
    grid_size: float = 1e-15,
):
    """Check if two polygons are not too close."""
    return poly1.distance(poly2) > distance + grid_size


def boundary_distance_req(
    poly1: shapely.Polygon | shapely.MultiPolygon,
    poly2: shapely.Polygon | shapely.MultiPolygon,
    distance: float,
    grid_size: float = 1e-15,
):
    """Check if the boundary of two polygons are not too close."""
    return poly1.boundary.distance(poly2.boundary) > distance + grid_size
