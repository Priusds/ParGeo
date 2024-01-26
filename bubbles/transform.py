"""Transform for the bubbles package."""
import math
from abc import ABC, abstractmethod
from typing import Any, Callable

from shapely import GeometryCollection, MultiPolygon, Polygon
from shapely import box as shapely_box
from shapely import union_all
from shapely.affinity import translate as shapely_translate

from .topology import Topology


class Transform(ABC):
    """Abstract class for transforms."""

    @abstractmethod
    def __call__(
        self, polygon: Polygon, level: int, topology: Topology
    ) -> Polygon | MultiPolygon:
        """Apply the transform to the polygon."""
        raise NotImplementedError


class Periodic(Transform):
    """Periodic transform."""

    def __init__(
        self,
        levels: (int | list[int]) | str = "any",
        x_length: float | list[float] = float("inf"),
        y_length: float | list[float] = float("inf"),
        alpha: float | list[float] = 0,
    ):
        """Define the Periodic transform.

        Args:
            levels: Levels to apply the periodicity. If "any" then the periodicity
                is applied to all levels. If a list of levels is given then the
                periodicity can be specified for each level individually.

            x_length: Periodic length in x direction. If "inf" then nothing
                happends in x direction. If a list is given then it must match the
                length of the levels list. Each entry is associated with the level
                with the corresponding index.

            y_length: Periodic length in y direction. Same as for
                `periodic_length_x`.

            alpha: Rotation angle in radians. Rotate the (x,y)-Plane ccw by alpha.
                NOTE: The polygons itself are not rotated, only the periodicity
                direction.
        """
        # Validate input
        if isinstance(levels, str) and not levels == "any":
            raise ValueError(
                f"Levels specification as string must be 'any' but given {levels}."
            )

        # Bring input in list form.
        if isinstance(levels, (str, int)):
            levels = [levels]
        if isinstance(x_length, (float, int)):
            x_length = [x_length]
        if isinstance(y_length, (float, int)):
            y_length = [y_length]
        if isinstance(alpha, (float, int)):
            alpha = [alpha]

        # Adjust lists
        if len(x_length) == 1:
            x_length = x_length * len(levels)
        if len(y_length) == 1:
            y_length = y_length * len(levels)
        if len(alpha) == 1:
            alpha = alpha * len(levels)

        # Check if all lists have the same length.
        if not len(set(map(len, [levels, x_length, y_length, alpha]))) == 1:
            raise ValueError("Input size mismatch.")

        # Check level uniqueness.
        if not len(set(levels)) == len(levels):
            raise ValueError("levels must not contain duplicates.")
        if "any" in levels and not len(levels) == 1:
            raise ValueError("If 'any' is given, no other levels can be given.")

        # Validate input types.
        for _input in zip(levels, x_length, y_length, alpha):
            self.__validate_input(*_input)

        self.any = "any" in levels

        if self.any:
            self.__lvl2xlength = {"any": x_length[0]}
            self.__lvl2ylength = {"any": y_length[0]}
            self.__lvl_2_alpha = {"any": alpha[0]}
        else:
            self.__lvl2xlength = dict(zip(levels, x_length))
            self.__lvl2ylength = dict(zip(levels, y_length))
            self.__lvl_2_alpha = dict(zip(levels, alpha))

    @property
    def affected_levels(self):
        """Return the levels affected by the periodicity transform."""
        return self.__lvl2xlength.keys()

    def __validate_input(
        self, level: int | str, x_length: float, y_length: float, alpha: float
    ):
        """Validate the input for the periodicity."""
        if not isinstance(level, int) and level != "any":
            raise ValueError(f"level must be int or `any` but given {level}.")
        if x_length <= 0:
            raise ValueError(f"x_length must be positive float but given {x_length}.")
        if y_length <= 0:
            raise ValueError(f"y_length must be positve float but given {y_length}.")
        if not isinstance(alpha, (float, int)):
            raise ValueError(f"alpha must be float but given {alpha}.")

    def update(self, level: int | str, x_length: float, y_length: float, alpha: float):
        """Update the periodic behavior of the transform.

        Note: If level is "any" then the periodicity is applied to all levels, meaning
        all levels are overwritten.
        """
        self.__validate_input(level, x_length, y_length, alpha)

        if self.any:
            self.__lvl2xlength = {"any": x_length[0]}
            self.__lvl2ylength = {"any": y_length[0]}
            self.__lvl_2_alpha = {"any": alpha[0]}
        else:
            self.__lvl2xlength[level] = x_length
            self.__lvl2ylength[level] = y_length
            self.__lvl_2_alpha[level] = alpha

    def __call__(
        self, polygon: Polygon, level: int, topology: Topology
    ) -> Polygon | MultiPolygon:
        """Apply the periodicity to the polygon."""
        # Check if the polygon is affected by the periodicity.
        if level not in self.affected_levels and not self.any:
            return polygon

        level = "any" if self.any else level

        x_length = self.__lvl2xlength[level]
        y_length = self.__lvl2ylength[level]
        alpha = self.__lvl_2_alpha[level]

        span = (
            topology.domain.bounds[2]
            - topology.domain.bounds[0]
            + topology.domain.bounds[3]
            - topology.domain.bounds[1]
        )
        x_reps = math.ceil(span / x_length)
        y_reps = math.ceil(span / y_length)

        x_dir = (
            [
                (math.cos(alpha), math.sin(alpha)),
                (-math.cos(alpha), -math.sin(alpha)),
            ]
            if x_reps > 0
            else []
        )
        y_dir = (
            [
                (math.sin(alpha), -math.cos(alpha)),
                (-math.sin(alpha), math.cos(alpha)),
            ]
            if y_reps > 0
            else []
        )

        # Pass the periodicity to the Repeat transform. As it is a special case of it.
        if x_reps > 0 and y_reps > 0:
            periodic_polygons = union_all(
                [
                    Repeat(v_dir, x_length, x_reps, w_dir, y_length, y_reps)(
                        polygon, level, topology
                    )
                    for v_dir in x_dir
                    for w_dir in y_dir
                ],
                grid_size=topology.grid_size,
            )
        elif x_reps > 0:
            periodic_polygons = union_all(
                [
                    Repeat(v_dir, x_length, x_reps, clip=False)(
                        polygon, level, topology
                    )
                    for v_dir in x_dir
                ],
                grid_size=topology.grid_size,
            )
        elif y_reps > 0:
            periodic_polygons = union_all(
                [
                    Repeat(w_dir, y_length, y_reps, clip=False)(
                        polygon, level, topology
                    )
                    for w_dir in y_dir
                ],
                grid_size=topology.grid_size,
            )
        else:
            return polygon

        return periodic_polygons


class Repeat(Transform):
    """Repeat transform."""

    def __init__(
        self,
        v_dir: tuple[float, float],
        v_scalar: float,
        v_rep: int,
        w_dir: tuple[float, float] = None,
        w_scalar: float = None,
        w_rep: int = None,
        clip: bool = True,
    ) -> Any:
        """Initialize a Repeat object.

        Args:
            v_dir: Direction vector. Will be normalized.
            v_scalar: Repeat distance.
            v_rep: Number of repeats.

            Note: Either all None or none of them.
            w_dir: Direction vector. Will be normalized.
            w_scalar: Repeat distance.
            w_rep: list[int]
        """
        if any((w_dir is None, w_scalar is None, w_rep is None)):
            if not all((w_dir is None, w_scalar is None, w_rep is None)):
                raise ValueError("Only")

        self.v_dir = v_dir
        self.v_scalar = v_scalar
        self.v_rep = v_rep

        self.w_dir = w_dir if w_dir is not None else (0, 0)
        self.w_scalar = w_scalar if w_scalar is not None else 0
        self.w_rep = w_rep if w_rep is not None else 0

        self.v_data = (v_dir, v_scalar, v_rep)
        self.w_data = (w_dir, w_scalar, w_rep) if w_dir is not None else None

        self.clip = clip

    def __call__(self, polygon: Polygon, level: int, topo: Topology):
        """Repeated the polygon in the defined directions.

        If `self.clip = True` then the repeated polygons will be clipped to the domain.
        """
        shifted_polygons = []

        bounds = polygon.bounds

        for i in range(self.v_rep + 1):
            for j in range(self.w_rep + 1):
                x_off = i * self.v_scalar * self.v_dir[0]
                y_off = i * self.v_scalar * self.v_dir[1]
                x_off += j * self.w_scalar * self.w_dir[0]
                y_off += j * self.w_scalar * self.w_dir[1]

                # If clip is set to True skip if the shifted polygon is
                # outside the domain.
                if self.clip:
                    bounds_shifted = (
                        bounds[0] + x_off,
                        bounds[1] + y_off,
                        bounds[2] + x_off,
                        bounds[3] + y_off,
                    )
                    # fast sufficient check if the shifted polygon is outside domain
                    if not shapely_box(*bounds_shifted).intersects(topo.domain):
                        continue

                shifted_polygons.append(
                    shapely_translate(polygon, xoff=x_off, yoff=y_off)
                )

        # Merge all shifted polygons.
        repeated_polygon = union_all(shifted_polygons)

        # Clip to domain.
        if self.clip:
            repeated_polygon = repeated_polygon.intersection(
                topo.domain, grid_size=topo.grid_size
            )
            if isinstance(repeated_polygon, GeometryCollection):
                repeated_polygon = MultiPolygon(
                    [poly for poly in repeated_polygon if poly.area > topo.grid_size]
                )

        return repeated_polygon


class Diffeomorphism(Transform):
    """Diffeomorphism transform."""

    def __init__(
        self,
        mapping: Callable[[list[float], list[float]], tuple[list[float], list[float]]],
    ):
        """Generate a diffeomorphism transform. The mapping is applyed to the
        boundary of the polygon.

        Args:
            mapping: x_coords, y_coords -> x_mapped, y_mapped
        """
        self.__mapping = mapping

    def __call__(
        self, polygon: Polygon | MultiPolygon, level: int, topology: Topology
    ) -> Polygon | MultiPolygon:
        """Apply the diffeomorphism to the polygon."""
        if isinstance(polygon, Polygon):
            return self.__map_polygon(polygon)
        elif isinstance(polygon, MultiPolygon):
            return MultiPolygon([self.__map_polygon(poly) for poly in polygon])

    def __map_polygon(self, polygon: Polygon):
        """Transform the polygon boundary."""
        x, y = polygon.exterior.xy
        x_mapped, y_mapped = self.__mapping(x, y)

        shell = [(xx, yy) for xx, yy in zip(x_mapped, y_mapped)]

        holes = []
        if len(polygon.interiors) > 0:
            for hole in polygon.interiors:
                x, y = hole.exterior
                x_mapped, y_mapped = self.__mapping(x, y)
                hole_mapped = [(xx, yy) for xx, yy in zip(x_mapped, y_mapped)]
                holes.append(hole_mapped)
        return Polygon(shell, holes)
