"""Utilities for type checking and type hinting."""
from typing import Union

from shapely.geometry import (
    GeometryCollection,
    LinearRing,
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)

# This is imported in the other modules.
GeometryCollection

SubDomain = Union[Polygon, MultiPolygon]
Level = int

Vector = tuple[float, float]


SHAPELY_GEOMETRIES = (
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
    LinearRing,
)

ShapelyGeometry = (
    Point
    | LineString
    | Polygon
    | MultiPoint
    | MultiLineString
    | MultiPolygon
    | GeometryCollection
    | LinearRing
)
