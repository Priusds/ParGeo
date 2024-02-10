"""
This is the main module of `pargeo`.

The module provides the `Topology` class, which can be used to generate complex
two-dimensional geometries, also called domains. 

The `Topology` class simplifies the process of sequentially modifing the domain
by adding geometries, in form of shapely polygons, to some initial domain.

The polygons that are added to the domain are equipped with a level, which is used 
to introduce a hierarchy between them. The higher the level of the polygon is,
the more visible it is, i.e. it covers the polygons with a lower level.

Additionally, holes can be defined, holes are levels, and polygons whose level
is in holes are not visible, basically they are holes in the domain.

Example:

>>> from pargeo.topology import Topology
>>> from shapely.geometry import Polygon
<BLANKLINE>
>>> domain = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
>>> topo = Topology(domain, holes={1})
>>> poly1 = Polygon(...) # Some shapely polygon
>>> topo.add(poly1, level=1)

"""

from __future__ import annotations

from collections import UserDict
from typing import Callable, Mapping, Optional, Sequence, Union

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from shapely import GeometryCollection, MultiPolygon, Polygon

DEFAULT_GRID_SIZE = 1e-15
INITIAL_LEVEL = 0


class Topology:
    """Class for the topology of the two-dimensional domain."""

    def __init__(
        self,
        domain: Polygon | MultiPolygon,
        holes: set[int] = set(),
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        """Initialize the topology.

        Args:
            domain: The domain.
            holes: Levels that are considered as holes.
                This can be changed at any time through `set_holes`.
            grid_size: The grid size.
        """
        if isinstance(domain, Polygon):
            domain = MultiPolygon([domain])
        self.__domain = domain
        self.__holes = holes  # Can be changed at any time.
        self.__grid_size = grid_size
        self.__lvl2multipoly = {INITIAL_LEVEL: domain}

    @property
    def holes(self):
        """Return the levels that are holes."""
        return self.__holes

    def set_holes(self, new_holes: set[int]):
        """Define the levels that are holes."""
        self.__holes = new_holes

    @property
    def grid_size(self) -> float:
        return self.__grid_size

    @property
    def domain(self) -> Polygon | MultiPolygon:
        return self.__domain

    @property
    def levels(self) -> list[int]:
        """List of sorted levels, without duplicates."""
        return sorted(self.__lvl2multipoly.keys())

    @property
    def lvl2multipoly(self):
        return self.__lvl2multipoly

    def add(
        self,
        polygon: Polygon | MultiPolygon,
        level: int,
        transform: Optional[
            Callable[
                [Polygon | MultiPolygon, int, Topology],
                MultiPolygon,
            ]
        ] = None,
        constraint: Optional[
            Callable[[Polygon | MultiPolygon, int, Topology], bool]
        ] = None,
        extend_domain: bool = False,
    ) -> bool:
        """
        Add a polygon to the topology.

        Args:
            polygon: The polygon to add.
            level: The level of the polygon.
            transform: Is called to transform the polygon before adding it.
                If will call `transform(polygon, level, self)`.
            constraint: Is called to check if the polygon satisfies the constraints.
                The constraints are checked after the polygon is transformed.
                If the constraints are not satisfied, the polygon is not added.
                It will call `constraints(polygon, level, self)`.
            extend_domain: Flag to extend the domain if the polygon is outside of it.
        """
        if level < 0:
            raise ValueError("`level` must be non-negative.")

        # Apply transform if given
        if transform is not None:
            polygon = transform(polygon, level, self)

        # Clip the transformed polygon w.r.t. the domain
        if not extend_domain:
            polygon = polygon.intersection(self.domain, grid_size=self.grid_size)
        else:
            extended_domain = self.domain.union(polygon, grid_size=self.grid_size)
            if isinstance(extended_domain, Polygon):
                extended_domain = MultiPolygon([extended_domain])
            self.__domain = extended_domain

        # Make sure polygon is a MultiPolygon with non-zero area polygons.
        if isinstance(polygon, (GeometryCollection, MultiPolygon)):
            polygon = MultiPolygon(
                [p for p in polygon.geoms if p.area > self.grid_size]
            )

        # Make sure polygon intersects with the reference domain
        if polygon.area < self.grid_size:
            return False

        # TODO: Apply segmentize strategy

        # Check if the polygon satisfies the constraints
        if constraint is not None:
            if not constraint(polygon, level, self):
                return False

        # Update self.__lvl2multipoly, this is the main logic of the class.
        # Multipolygons with a higher level are cutted out from those with a lower
        # level, if the intersect. The level is a visibility ranking.
        new_lvl2multipoly = dict()
        for running_level, running_polygon in self.__lvl2multipoly.items():
            if not polygon.intersects(running_polygon):
                # running_polygon does not intersect with `polygon` so it can be
                # added directly to the new `__lvl2multipoly`.
                new_lvl2multipoly[running_level] = running_polygon
                continue

            # running_polygon intersects with `polygon`. Both might
            # need to be updated.
            if level != running_level:
                higher_polygon = polygon if level > running_level else running_polygon
                lower_polygon = polygon if level < running_level else running_polygon

                # Cut out the higher level polygon from the lower level polygon
                lower_polygon = lower_polygon.difference(
                    higher_polygon, grid_size=self.grid_size
                )

                # After shapely.difference, lower_polygon can be a GeometryCollection.
                if isinstance(lower_polygon, GeometryCollection):
                    lower_polygon = MultiPolygon(
                        [p for p in lower_polygon.geoms if p.area > self.grid_size]
                    )

                # If the lower polygon is approximately empty, it is not added.
                # If polygon = lower_polygon, then the polygon is not added.
                if lower_polygon.area < self.grid_size:
                    if level < running_level:
                        return False
                    continue

                # Add the intersection points, this is needed to create delauny
                # triangulations by GMSH.
                intersection_ = higher_polygon.intersection(
                    lower_polygon, grid_size=self.grid_size
                )

                # `intersection` can be a GeometryCollection, if so, iterate over the
                # geoms and add them one by one to the higher_polygon.
                if isinstance(intersection_, GeometryCollection):
                    for geom in intersection_.geoms:
                        # TODO: Check why geom can be a polygon, sometimes
                        # this happens, and the polygon has some approximate 0 area.

                        # Only add geoms with 0 area.
                        if geom.area == 0:
                            # Note geom.area == 0 is needed to remove polygons that,
                            # given our grid_size are not visible, *but* if a new
                            # polygon is added with some epsilon area, then this needs
                            # to be removed, which happens in the next if statement.
                            higher_polygon = higher_polygon.union(
                                geom, grid_size=self.grid_size
                            )
                            if isinstance(higher_polygon, GeometryCollection):
                                higher_polygon = MultiPolygon(
                                    [
                                        poly
                                        for poly in higher_polygon.geoms
                                        if poly.area > self.grid_size
                                    ]
                                )
                else:
                    higher_polygon = higher_polygon.union(
                        intersection_, grid_size=self.grid_size
                    )
                # TODO: Is it possible that the union in the upper if statement
                # is again a Polygon? If so, then we need to convert it.
                if isinstance(higher_polygon, Polygon):
                    higher_polygon = MultiPolygon([higher_polygon])

                # TODO: Is this necessary?
                if isinstance(higher_polygon, GeometryCollection):
                    higher_polygon = MultiPolygon(
                        [p for p in higher_polygon.geoms if p.area > self.grid_size]
                    )

                # Update the running level polygon and polygon.
                if running_level < level:
                    assert not isinstance(lower_polygon, GeometryCollection)
                    new_lvl2multipoly[running_level] = lower_polygon
                    polygon = higher_polygon  # polygon gets bigger

                else:  # running_level > level
                    assert not isinstance(higher_polygon, GeometryCollection)
                    new_lvl2multipoly[running_level] = higher_polygon
                    polygon = lower_polygon  # polygon gets smaller

        # Merge the polygon with the multipolygon of the same level.
        if level in self.__lvl2multipoly:
            polygon = polygon.union(
                self.__lvl2multipoly[level], grid_size=self.grid_size
            )
            # The shapely.union might create a GeometryCollection.
            # Appearently, polygon.union(multipolygon) can create
            # a GeometryCollection with a multiLineString.
            if isinstance(polygon, GeometryCollection):
                polygon = MultiPolygon(
                    [p for p in polygon.geoms if p.area > self.grid_size]
                )

        if isinstance(polygon, Polygon):
            polygon = MultiPolygon([polygon])

        if not isinstance(polygon, MultiPolygon):
            raise ValueError(f"`polygon` is not a MultiPolygon, but a {type(polygon)}.")

        new_lvl2multipoly[level] = polygon
        self.__lvl2multipoly = new_lvl2multipoly
        return True

    def plot(self, color_holes: bool = False, color_hole_boundaries=True) -> None:
        """Plot the topology."""
        colormap = DefaultColors.get_color_map(self.levels, self.holes, color_holes)

        def plot_tree_rec(node: Node) -> None:
            if not node.is_root:
                # plot node
                x, y = node["polygon"].exterior.xy
                plt.fill(x, y, color=colormap[node["level"]])

                if node["level"] not in self.holes:
                    plt.plot(x, y, "-", linewidth=1, color=DefaultColors.boundary)
                if color_hole_boundaries and node["level"] in self.holes:
                    plt.plot(x, y, "--", linewidth=0.5, color=DefaultColors.boundary)

                # plot interiors in white
                for interior in node["polygon"].interiors:
                    x, y = interior.xy
                    plt.fill(x, y, color="white")
                    plt.plot(x, y, "-", linewidth=1, color="black")

            for child in node["children"]:
                plot_tree_rec(child)

        make_legend(self.holes, colormap)

        plot_tree_rec(self.get_inclusion_tree().root)
        plt.show()

    def flatten(self) -> list[tuple[Polygon, int]]:
        """Return a list of polygons and their levels."""
        flattened_polygons = []
        for level, running_polygon in self.__lvl2multipoly.items():
            if isinstance(running_polygon, MultiPolygon):
                for polygon in running_polygon.geoms:
                    flattened_polygons.append((polygon, level))
            else:
                flattened_polygons.append((running_polygon, level))
        return flattened_polygons

    def get_inclusion_tree(self) -> InclusionTree:
        """Return the topology as a tree."""
        root_node = Node(level=INITIAL_LEVEL - 1, polygon=Polygon(), children=[])
        tree = InclusionTree(root=root_node)

        for polygon, level in self.flatten():
            tree.add(polygon, level)
        return tree


class InclusionTree:
    """Represent the flattened topology as a tree.

    The tree is used to plot the topology.

    The tree is a directed acyclic graph, where the root is just an auxiliary
    node. Each (non-root) node represents a polygon, with a level. All polygons
    represented by the nodes are disjoint. The nodes represent the polygons in
    the flattened topology.

    This tree encodes the semi-ordering of the polygons, that is given by the
    inclusion relation. This means, a directed edge from node A to node B exists
    if and only if the polygon represented by node A includes the polygon of node B.
    Note: The polygon of the exterior of the polygon (spoken in terms of shapely) is
    taken into account for the inclusion.
    """

    def __init__(self, root: Node) -> None:
        root.is_root = True
        self.__root = root

    @property
    def root(self) -> Node:
        return self.__root

    def add(self, polygon: Polygon, level: int) -> None:
        """
        Recursively add a node to the tree.

        This function adds a node to the tree by finding its correct position in relation to the parent node.
        The node is added as a child of the parent node if the parent node includes the node.
        If the node includes a child of the parent node, the child is removed from the parent node's children
        and added to the node's children. The depth of the node and any affected children is updated accordingly.

        Args:
            parent: The parent node.
            node: The node to be added.

        Returns:
            None
        """

        def add_node_rec(parent: Node, node: Node) -> None:
            for child in parent["children"]:
                if child.includes(node):
                    node.depth = child.depth + 1
                    add_node_rec(child, node)
                    return

                elif node.includes(child):
                    node.depth = child.depth
                    child.depth = node.depth + 1
                    node["children"].append(child)
                    parent["children"].remove(child)
                    add_node_rec(parent, node)
                    return
                else:
                    continue
            node.depth = parent.depth + 1
            parent["children"].append(node)

        node = Node(level=level, polygon=polygon)
        add_node_rec(self.__root, node)

    def show(self) -> None:
        """Print the tree recursively."""

        def show_rec(node: Node, depth: int) -> None:
            if not node.is_root:
                print("    " * (depth - 1), node)
            for child in node["children"]:
                show_rec(child, depth + 1)

        show_rec(self.__root, 0)


class Node(UserDict):
    """Node of the inclusion tree, representing a polygon.

    Nodes are semi-ordered by the inclusion relation.
    """

    def __init__(self, level: int, polygon: Polygon, children=None):
        super().__init__()
        self.data = {
            "level": level,
            "polygon": polygon,
            "children": children if children is not None else [],
        }
        self.is_root = False
        self.depth = 0

    def includes(self, node: Node) -> bool:
        """Check if this node includes the incoming node."""
        return Polygon(self["polygon"].exterior).contains(node["polygon"])

    def __repr__(self) -> str:
        return f"Node(level={self['level']}, depth={self.depth}), n_children={len(self['children'])}"


Color = Union[str, tuple[int, int, int], tuple[int, int, int, int]]


class DefaultColors:
    """Default colors for plotting."""

    boundary = "black"
    domain_boundary = "black"
    domain = "silver"
    hole = "white"
    interior_filling = "white"
    inter_hole_bound = "white"

    @staticmethod
    def get_color_map(
        levels: Sequence[int], holes: set[int], color_holes=False
    ) -> Mapping[int, Color]:
        """Get a color map for the levels."""
        lvl2cl = dict()
        level_set = set(levels)
        if 0 in level_set:
            lvl2cl[0] = DefaultColors.domain
            level_set.remove(0)
        if len(level_set) > 0:
            norm = Normalize(vmin=min(level_set), vmax=max(level_set))
            for lvl in level_set:
                if lvl in holes and not color_holes:
                    lvl2cl[lvl] = DefaultColors.hole
                else:
                    lvl2cl[lvl] = plt.cm.cool(norm(lvl))  # type: ignore
        return lvl2cl


def make_legend(holes: Sequence[int], colormap: Mapping[int, Color]):
    """Make a legend."""
    handles = []
    descriptions = []
    for lvl, color in colormap.items():
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=color, edgecolor="black"))
        if lvl in holes:
            descriptions.append(f"{lvl}. level (hole)")
        else:
            descriptions.append(f"{lvl}. level")
    plt.legend(
        handles,
        descriptions,
        loc="upper left",
        bbox_to_anchor=(0.95, 1),
    ).set_draggable(True)
