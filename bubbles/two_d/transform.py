

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