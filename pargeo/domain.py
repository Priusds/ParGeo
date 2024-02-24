"""
This is the main module of `pargeo`.

The module provides the `Domain` class, which can be used to generate complex
two-dimensional geometries, also called domains. 

The `Domain` class simplifies the process of sequentially modifing the domain
by adding geometries, in form of shapely polygons, to some initial domain.

The polygons that are added to the domain are equipped with a level, which is used 
to introduce a hierarchy between them. The higher the level of the polygon is,
the more visible it is, i.e. it covers the polygons with a lower level.

Additionally, holes can be defined, holes are levels, and polygons whose level
is in holes are not visible, basically they are holes in the domain.

Example:

>>> from pargeo.domain import Domain
>>> from shapely.geometry import Polygon
<BLANKLINE>
>>> domain = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
>>> domain = Domain(domain, holes={1})
>>> subdomain = Polygon(...) # Some shapely polygon
>>> domain.add(subdomain, level=1)

"""

from __future__ import annotations

from collections import UserDict
from typing import Any, Mapping, Protocol

import matplotlib.pyplot as plt

from shapely import GeometryCollection

from pargeo.utils.constants import BACKGROUND_LEVEL, DEFAULT_GRID_SIZE
from pargeo.utils.plot_utils import plot_legend, plot_polygon, DefaultColors, Color
from pargeo.utils.typing_utils import Level, MultiPolygon, Polygon, SubDomain


class Transform(Protocol):
    """
    The Transform protocol represents a callable that takes a SubDomain, a Level,
    a Domain, and any number of additional keyword arguments, and returns a SubDomain.

    This is used to define the interface for transformations that can be applied to a SubDomain.
    """

    def __call__(
        self, subdomain: SubDomain, level: Level, domain: Domain, **kwargs
    ) -> SubDomain:
        """
        Applies the transformation to the given SubDomain.

        Args:
            subdomain (SubDomain): The SubDomain to transform.
            level (Level): The Level of the subdomain.
            domain (Domain): The underlying domain.
            **kwargs: Additional keyword arguments for the transformation.

        Returns:
            SubDomain: The transformed SubDomain.
        """
        ...


class Constraint(Protocol):
    """
    The Constraint protocol represents a callable that takes a SubDomain, a Level,
    a Domain, and any number of additional keyword arguments, and returns a boolean.

    This is used to define the interface for constraints that can be applied to a SubDomain.
    """

    def __call__(
        self, subdomain: SubDomain, level: Level, domain: Domain, **kwargs
    ) -> bool:
        """
        Applies the constraint to the given SubDomain.

        Args:
            subdomain (SubDomain): The SubDomain to apply the constraint to.
            level (Level): The Level of the SubDomain.
            domain (Domain): The underlying domain.
            **kwargs: Additional keyword arguments for the constraint.

        Returns:
            bool: True if the SubDomain satisfies the constraint, False otherwise.
        """
        ...


class Domain:
    """Class for the creation of complex two-dimensional domains.

    A domain is a two-dimensional geometry, which can be modified by adding
    subdomains to it. The subdomains are shapely polygons, or multipolygons,
    and can be added using the `add` method. The subdomains are equipped with
    a level, which is used to introduce a hierarchy between them.

    Subdomains with a higher level will be cut out from those with a lower level.
    Visually speaking the higher level subdomains are more visible and cover the
    lower level subdomains. The lowest level is the `Domain.background_level`.

    Subdomains with the same level are merged together.
    """

    # Level of the background.
    background_level = BACKGROUND_LEVEL

    def __init__(
        self,
        background: SubDomain,
        holes: set[Level] = set(),
        grid_size: float = DEFAULT_GRID_SIZE,
    ) -> None:
        """Initialize the domain.

        Args:
            background (SubDomain): This subdomain is the background of the domain
                and is associated with the level `Domain.background_level`.

            holes (set[Level], optional): Set of levels that are considered as holes.
                Can be changed at any time using the `set_holes` method. Defaults to `set()`.

            grid_size (float, optional): Used `grid_size` for shapely geometry manipulating operations.
        """
        # Make sure background is a SubDomain
        if not isinstance(background, (Polygon, MultiPolygon)):
            raise ValueError(
                f"`background` must be a shapely Polygon or MultiPolygon, but is {type(background)}."
            )
        # Internally only work with MultiPolygons
        if isinstance(background, Polygon):
            background = MultiPolygon([background])
        self.__profile = background
        self.__holes = holes  # Can be changed at any time.
        self.__grid_size = grid_size
        self.__level_to_subdomain = {Domain.background_level: background}

    @property
    def level_to_subdomain(self) -> dict[Level, SubDomain]:
        """Dictionary mapping each level to its corresponding subdomain."""
        return self.__level_to_subdomain

    @property
    def holes(self) -> set[Level]:
        """Set of levels that are considered as holes."""
        return self.__holes

    def set_holes(self, new_holes: set[Level]):
        """Set `self.holes`."""
        self.__holes = new_holes

    @property
    def grid_size(self) -> float:
        """Used `grid_size` for shapely geometry manipulating operations."""
        return self.__grid_size

    @property
    def profile(self) -> SubDomain:
        """The profile of the domain.

        If the added subdomains are clipped to the background, then the profile
        matches the background. Otherwise it can be bigger.
        """
        return self.__profile

    @property
    def levels(self) -> list[Level]:
        """Present levels in sorted order."""
        return sorted(self.__level_to_subdomain.keys())

    @property
    def subdomains(self) -> list[SubDomain]:
        """Return a list of all subdomains."""
        subdomains = [subdomain for subdomain, _ in self.as_list()]
        return subdomains

    def add_subdomain(
        self,
        subdomain: SubDomain,
        level: Level,
        transform: Transform | None = None,
        constraint: Constraint | None = None,
        clip: bool = True,
        **kwargs: Any,
    ) -> bool:
        """
        Adds a subdomain to the domain.

        This method performs the following steps before potentially adding the subdomain:
            1. Transforms the incoming subdomain using `transform`.
            2. If clip is `True`, clips the subdomain to the domain's profile.
            3. If `constraint` are given, only proceede if these are met.
            4. Add the subdomain.

        Args:
            subdomain (SubDomain): The subdomain to be added.
            level (Level): The associated level of the subdomain.
            transform (Transform): A function to transform the subdomain before adding it.
                It will call `transform(subdomain, level, self)`.
            constraint (Constraint): A function to check if the subdomain should be added.
                It will call `constraints(subdomain, level, self)`.
            clip (bool): If `True`, clips the subdomain to the background.

        Note: If `level == Domain.background_level` and `clip == True`, then the
            background might be modified.

        Returns:
            bool: True if the subdomain was added, False otherwise.
        """
        if level < Domain.background_level:
            raise ValueError(
                f"`level` must be greater-equal to {Domain.background_level}."
            )

        # Apply transform if given
        if transform is not None:
            transform_kwargs = self._extract_kwargs(kwargs, "transform_")
            transform_kwargs["clip"] = clip
            subdomain = transform(subdomain, level, self, **transform_kwargs)

        # Clip the transformed subdomain to the background
        if clip:
            subdomain = subdomain.intersection(self.profile, grid_size=self.grid_size)
        else:
            profile = self.profile.union(subdomain, grid_size=self.grid_size)
            if isinstance(profile, Polygon):
                profile = MultiPolygon([profile])
            self.__profile = profile

        # Make sure subdomain is a MultiPolygon with positive area polygons
        if isinstance(subdomain, (GeometryCollection, MultiPolygon)):
            subdomain = MultiPolygon(
                [p for p in subdomain.geoms if p.area > self.grid_size]
            )

        # Make sure subdomain is not empty
        if subdomain.area < self.grid_size:
            return False

        # TODO: Apply segmentize strategy

        # Check constraints if given
        if constraint is not None:
            constraint_kwargs = self._extract_kwargs(kwargs, "constraint_")
            if not constraint(subdomain, level, self, **constraint_kwargs):
                return False

        # Update `self.__level_to_subdomain`
        level_to_subdomain = dict()
        for running_level, running_subdomain in self.__level_to_subdomain.items():
            if not subdomain.intersects(running_subdomain):
                # running_subdomain does not intersect with `subdomain` so it can be
                # added directly to `level_to_subdomain`.
                level_to_subdomain[running_level] = running_subdomain
                continue

            # `running_subdomain` intersects with `subdomain`
            # Both might need to be updated
            if level != running_level:
                higher_subdomain = (
                    subdomain if level > running_level else running_subdomain
                )
                lower_subdomain = (
                    subdomain if level < running_level else running_subdomain
                )

                # Cut out the higher level subdomain from the lower level subdomain
                lower_subdomain = lower_subdomain.difference(
                    higher_subdomain, grid_size=self.grid_size
                )

                # After shapely.difference, lower_subdomain can be a GeometryCollection.
                if isinstance(lower_subdomain, GeometryCollection):
                    lower_subdomain = MultiPolygon(
                        [p for p in lower_subdomain.geoms if p.area > self.grid_size]
                    )

                # If the lower subdomain is approximately empty, it is not added.
                # If subdomain = lower_subdomain, then the subdomain is not added.
                if lower_subdomain.area < self.grid_size:
                    if level < running_level:
                        return False
                    continue

                # Add the intersection points, this is needed to create delauny
                # triangulations by GMSH.
                intersection_ = higher_subdomain.intersection(
                    lower_subdomain, grid_size=self.grid_size
                )

                # `intersection` can be a GeometryCollection, if so, iterate over the
                # geoms and add them one by one to the higher_subdomain.
                if isinstance(intersection_, GeometryCollection):
                    for geom in intersection_.geoms:
                        # TODO: Check why geom can be a `shapely.geometry.Polygon`.
                        # This happens sometimes, and the `Polygon` has some approximate area 0.

                        # Only add shapely geometries with area == 0, e.g. Points
                        if geom.area == 0:
                            # Note `geom.area == 0` is needed to remove `Polygons` that,
                            # w.r.t. the given `grid_size` are not visible, *but* if a new
                            # `Polygon` is added with some epsilon area, then this needs
                            # to be removed, which happens in the next if statement.
                            higher_subdomain = higher_subdomain.union(
                                geom, grid_size=self.grid_size
                            )
                            if isinstance(higher_subdomain, GeometryCollection):
                                higher_subdomain = MultiPolygon(
                                    [
                                        poly
                                        for poly in higher_subdomain.geoms
                                        if poly.area > self.grid_size
                                    ]
                                )
                else:
                    higher_subdomain = higher_subdomain.union(
                        intersection_, grid_size=self.grid_size
                    )
                # TODO: Is it possible that the union in the upper if statement
                # is again a `Polygon`? If so, then we need to convert it.
                if isinstance(higher_subdomain, Polygon):
                    higher_subdomain = MultiPolygon([higher_subdomain])

                # TODO: Is this necessary?
                if isinstance(higher_subdomain, GeometryCollection):
                    higher_subdomain = MultiPolygon(
                        [p for p in higher_subdomain.geoms if p.area > self.grid_size]
                    )

                # Update  `subdomain` and `level_to_subdomain`
                if running_level < level:
                    assert not isinstance(lower_subdomain, GeometryCollection)
                    level_to_subdomain[running_level] = lower_subdomain
                    subdomain = higher_subdomain  # subdomain grows

                else:  # running_level > level
                    assert not isinstance(higher_subdomain, GeometryCollection)
                    level_to_subdomain[running_level] = higher_subdomain
                    subdomain = lower_subdomain  # subdomain shrinks

        # Merge the subdomain with the multipolygon of the same level.
        if level in self.__level_to_subdomain:
            subdomain = subdomain.union(
                self.__level_to_subdomain[level], grid_size=self.grid_size
            )
            # The shapely.union might create a GeometryCollection.
            # Appearently, polygon.union(multipolygon) can create
            # a GeometryCollection with a multiLineString.
            if isinstance(subdomain, GeometryCollection):
                subdomain = MultiPolygon(
                    [p for p in subdomain.geoms if p.area > self.grid_size]
                )

        if isinstance(subdomain, Polygon):
            subdomain = MultiPolygon([subdomain])

        if not isinstance(subdomain, MultiPolygon):
            raise ValueError(
                f"`subdomain` is not a MultiPolygon, but a {type(subdomain)}."
            )

        level_to_subdomain[level] = subdomain
        self.__level_to_subdomain = level_to_subdomain
        return True

    def plot(
        self,
        title: str = "Domain",
        color_holes: bool = False,
        color_hole_boundaries: bool = True,
        make_legend: bool = True,
        safe_file: str | None = None,
        color_map: Mapping[Level, Color] | None = None,
    ) -> None:
        """Plot the domain.

        Args:
            title (str, optional): Title of the plot. Defaults to "Domain".
            color_holes (bool, optional): If `True`, the holes are colored. Defaults to `False`.
            color_hole_boundaries (bool, optional): If `True`, the boundaries of the holes are colored.
                Defaults to `True`.
            make_legend (bool, optional): If `True`, a legend is added to the plot. Defaults to `True`.
            safe_file (str, optional): If given, the plot is saved to the file. Defaults to `None`.
            color_map (Mapping[Level, Color], optional): A mapping of levels to colors.
                If `None`, the default color map is used. Defaults to `None`.
        """
        if color_map is None:
            colormap = DefaultColors.get_color_map(
                self.levels, self.holes, color_holes, self.background_level
            )

        def plot_tree_rec(node: Node) -> None:
            """Recursively plot the domain."""
            if not node.is_root:
                node.plot(colormap, color_hole_boundaries, self.holes)

            for child in node["children"]:
                plot_tree_rec(child)

        plt.title(title)
        if make_legend:
            plot_legend(self.holes, colormap, Domain.background_level)
        plot_tree_rec(self.as_tree().root)
        if safe_file is not None:
            plt.savefig(safe_file, dpi=300)
        plt.show()

    def as_list(self) -> list[tuple[Polygon, Level]]:
        """Return a list of all subdomains with their levels."""
        polyongs_levels = []
        for level, running_polygon in self.__level_to_subdomain.items():
            if isinstance(running_polygon, MultiPolygon):
                for polygon in running_polygon.geoms:
                    polyongs_levels.append((polygon, level))
            else:
                polyongs_levels.append((running_polygon, level))
        return polyongs_levels

    def as_tree(self) -> InclusionTree:
        """Represent the domain as a tree.

        Each node of the tree represents a subdomain of the domain. If node A
        is a parent of node B, then the subdomain of node A includes the subdomain
        of node B.
        """
        root_node = Node(
            level=Domain.background_level - 1, polygon=Polygon(), children=[]
        )
        tree = InclusionTree(root=root_node)

        for polygon, level in self.as_list():
            tree.add(polygon, level)
        return tree

    def _extract_kwargs(self, kwargs: dict[str, Any], prefix: str) -> dict[str, Any]:
        """Extract the kwargs that start with `prefix`."""
        _kwargs_id = prefix
        _n = len(_kwargs_id)
        return {k[_n]: v for k, v in kwargs.items() if k.startswith(_kwargs_id)}


class InclusionTree:
    """Represents a `Domain` object as a tree.

    The tree is used to plot the domain.

    The tree is a directed acyclic graph, where the root is just an auxiliary
    node. Each (non-root) node represents a polygon, with a level. All polygons
    represented by the nodes are disjoint. The nodes represent the polygons in
    the flattened domain.

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

    def add(self, polygon: Polygon, level: Level) -> None:
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
            """Recursively add a node to the tree."""
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
            """Recursively print each node."""
            if not node.is_root:
                print("    " * (depth - 1), node)
            for child in node["children"]:
                show_rec(child, depth + 1)

        show_rec(self.__root, depth=0)


class Node(UserDict):
    """Node of the inclusion tree, representing a polygon.

    Nodes are semi-ordered by the inclusion relation.
    """

    def __init__(self, level: Level, polygon: Polygon, children=None):
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

    def plot(
        self,
        colormap: Mapping[Level, Color],
        color_hole_boundaries: bool,
        holes: set[int],
        show: bool = False,
        boundary_color: Color | None = None,
    ):
        """Plot the node."""
        is_hole = self.data["level"] in holes
        color = colormap[self.data["level"]]
        polygon = self.data["polygon"]
        if boundary_color is None:
            boundary_color = DefaultColors.boundary

        if not is_hole:
            plot_polygon(polygon, color, boundary_color, "-", 1.0, show)
        else:
            if color_hole_boundaries:
                plot_polygon(polygon, color, boundary_color, "--", 0.5, show)
            else:
                plot_polygon(polygon, color, "none", "-", 1.0, show)
