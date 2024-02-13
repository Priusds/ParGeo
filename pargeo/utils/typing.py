from typing import Union

from shapely.geometry import (GeometryCollection, LinearRing, LineString,
                              MultiLineString, MultiPoint, MultiPolygon, Point,
                              Polygon)

# This is imported in the other modules.
GeometryCollection

SubDomain = Union[Polygon, MultiPolygon]
Level = int

Vector = tuple[float, float]
Color = Union[str, tuple[int, int, int], tuple[int, int, int, int]]

Shapely_Geometry = Union[
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
    LinearRing,
]
