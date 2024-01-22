import shapely
from shapely.affinity import translate as shapely_translate

class Transform(object):
    def __call__(
        self, polygon: shapely.Polygon | shapely.MultiPolygon
    ) -> shapely.Polygon | shapely.MultiPolygon:
        pass


class PeriodicXY(Transform):
    def __init__(self, midpoint, width, height):
        pass


class PeriodicX(Transform):
    def __init__(self, midpoint, width, height):
        pass


class PeriodicY(Transform):
    def __init__(self, midpoint, width, height):
        pass


def clip_x(
    polygon: shapely.Polygon | shapely.MultiPolygon,
    x_min: float,
    x_max: float,
    grid_size: float = 1e-15,
):
    """Clip the polygon along the x-axis."""
    poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
    if poly_x_min < x_min - grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(x_min, poly_y_min, poly_x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, poly_y_min, x_min, poly_y_max),
                grid_size=grid_size,
            ),
            xoff=x_max - x_min,
        )
    elif poly_x_max > x_max + grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, poly_y_min, x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(x_max, poly_y_min, poly_x_max, poly_y_max),
                grid_size=grid_size,
            ),
            xoff=x_min - x_max,
        )
    else:
        assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
        return polygon
    if isinstance(poly_translated, shapely.GeometryCollection):
        poly_translated = shapely.MultiPolygon(
            [p for p in poly_translated.geoms if p.area > grid_size]
        )
    polygon = poly_in.union(poly_translated, grid_size=grid_size)
    assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
    return clip_x(polygon, x_min, x_max, grid_size=grid_size)


def clip_y(
    polygon: shapely.Polygon | shapely.MultiPolygon,
    y_min: float,
    y_max: float,
    grid_size: float = 1e-15,
):
    """Clip the polygon along the y-axis."""
    poly_x_min, poly_y_min, poly_x_max, poly_y_max = polygon.bounds
    if poly_y_min < y_min - grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, y_min, poly_x_max, poly_y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, poly_y_min, poly_y_max, y_min),
                grid_size=grid_size,
            ),
            yoff=y_max - y_min,
        )
    elif poly_y_max > y_max + grid_size:
        poly_in = shapely.intersection(
            polygon,
            shapely.geometry.box(poly_x_min, poly_y_min, poly_x_max, y_max),
            grid_size=grid_size,
        )
        poly_translated = shapely_translate(
            shapely.intersection(
                polygon,
                shapely.geometry.box(poly_x_min, y_max, poly_x_max, poly_y_max),
                grid_size=grid_size,
            ),
            yoff=y_min - y_max,
        )
    else:
        assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
        return polygon
    if isinstance(poly_translated, shapely.GeometryCollection):
        poly_translated = shapely.MultiPolygon(
            [p for p in poly_translated.geoms if p.area > grid_size]
        )
    polygon = poly_in.union(poly_translated, grid_size=grid_size)
    assert isinstance(polygon, (shapely.Polygon, shapely.MultiPolygon))
    return clip_y(polygon, y_min, y_max, grid_size=grid_size)
