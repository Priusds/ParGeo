"""Transform for the pargeo package."""
import math
from typing import Any, Callable

from pargeo.domain import Domain, Transform
from pargeo.utils.geometry_utils import repeat
from pargeo.utils.typing_utils import Level, MultiPolygon, Polygon, SubDomain, Vector


class Periodic(Transform):
    """Periodic transform."""

    def __init__(
        self,
        levels: int | list[int] | str = "all",
        x_length: float | list[float] = float("inf"),
        y_length: float | list[float] = float("inf"),
        alpha: float | list[float] = 0,
    ):
        """Define the Periodic transform.

        Each level can have its own periodicity. If levels is "all" then all levels
        are treaded equally and the periodicity is applied to all levels.

        Args:
            levels: Levels to apply the periodicity. If "all" then the periodicity
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
        levels_ = levels if isinstance(levels, list) else []
        x_length = x_length if isinstance(x_length, list) else [x_length]
        y_length = y_length if isinstance(y_length, list) else [y_length]
        alpha = alpha if isinstance(alpha, list) else [alpha]
        self.all = False
        self.all_xlength = 0.0
        self.all_ylength = 0.0
        self.all_alpha = 0.0

        self.__lvl2xlength = {}
        self.__lvl2ylength = {}
        self.__lvl_2_alpha = {}

        # Validate input
        n_levels = 0
        if isinstance(levels, str):
            if levels == "all":
                if not (len(x_length) == 1 and len(y_length) == 1 and len(alpha) == 1):
                    raise ValueError(
                        "If levels is 'all' then x_length, y_length and alpha must be"
                        "floats."
                    )
                self.all = True
                n_levels = 1
                self.all_xlength = x_length[0]
                self.all_ylength = y_length[0]
                self.all_alpha = alpha[0]
            else:
                raise ValueError(
                    f"Levels specification as string must be 'all' but given {levels}."
                )
        elif isinstance(levels, int):
            levels_.append(levels)
            if not (len(x_length) == 1 and len(y_length) == 1 and len(alpha) == 1):
                raise ValueError(
                    "If levels is 'all' then x_length, y_length and alpha must be"
                    "floats."
                )
            self.all = False
            n_levels = 1
        elif isinstance(levels, list):
            n_levels = len(levels)

            x_length = x_length * n_levels if len(x_length) == 1 else x_length
            y_length = y_length * n_levels if len(y_length) == 1 else y_length
            alpha = alpha * n_levels if len(alpha) == 1 else alpha
            self.all = False
        else:
            raise ValueError(
                f"levels must be int, list[int] or 'all' but given {levels}."
            )
        if not self.all:
            self.__lvl2xlength = dict(zip(levels_, x_length))
            self.__lvl2ylength = dict(zip(levels_, y_length))
            self.__lvl_2_alpha = dict(zip(levels_, alpha))

    @property
    def affected_levels(self):
        """Return the levels affected by the periodicity transform."""
        return self.__lvl2xlength.keys()

    def __validate_input(
        self, level: int | str, x_length: float, y_length: float, alpha: float
    ):
        """Validate the input for the periodicity."""
        if not isinstance(level, int) and level != "all":
            raise ValueError(f"level must be int or `all` but given {level}.")
        if x_length <= 0:
            raise ValueError(f"x_length must be positive float but given {x_length}.")
        if y_length <= 0:
            raise ValueError(f"y_length must be positve float but given {y_length}.")
        if not isinstance(alpha, (float, int)):
            raise ValueError(f"alpha must be float but given {alpha}.")

    def update(self, level: int, x_length: float, y_length: float, alpha: float):
        """Update the periodic behavior of the transform.

        Note: If level is "all" then the periodicity is applied to all levels, meaning
        all levels are overwritten.
        """
        self.__validate_input(level, x_length, y_length, alpha)
        self.__lvl2xlength[level] = x_length
        self.__lvl2ylength[level] = y_length
        self.__lvl_2_alpha[level] = alpha

    def __call__(
        self, subdomain: SubDomain, level: int, domain: Domain, **kwargs: Any
    ) -> SubDomain:
        """Apply the periodicity to the polygon."""
        # Check if the polygon is affected by the periodicity.
        if level not in self.affected_levels and not self.all:
            return subdomain

        x_length = self.__lvl2xlength[level] if not self.all else self.all_xlength
        y_length = self.__lvl2ylength[level] if not self.all else self.all_ylength
        alpha = self.__lvl_2_alpha[level] if not self.all else self.all_alpha

        span: float = (
            domain.profile.bounds[2]
            - domain.profile.bounds[0]
            + domain.profile.bounds[3]
            - domain.profile.bounds[1]
        )
        n_reps_x = math.ceil(span / x_length)
        n_reps_y = math.ceil(span / y_length)

        return repeat(
            subdomain,
            ((math.cos(alpha), math.sin(alpha)), x_length, n_reps_x, True),
            ((math.sin(alpha), -math.cos(alpha)), y_length, n_reps_y, True),
            clipping_plane=domain.profile,
            grid_size=domain.grid_size,
        )


class Repeat(Transform):
    """Repeat transform."""

    def __init__(
        self,
        repeat_info: tuple[Vector, float, int, bool],
        repeat_info_2: tuple[Vector, float, int, bool] | None = None,
        excluded_levels: list[Level] = [],
    ) -> None:
        """Initialize a Repeat object.

        Args:
            repeat_info: Some info
            repeat_info_2: Some info
            excluded_levels: Levels to exclude from the transformation.
        """
        if repeat_info_2 is None:
            repeat_info_2 = ((0, 0), 0, 0, False)

        self.repeat_info = repeat_info
        self.repeat_info_2 = repeat_info_2
        self.excluded_levels = excluded_levels

    def __call__(self, subdomain: SubDomain, level: int, domain: Domain, **kwargs: Any):
        """Apply the repeat transform to the SubDomain.

        Args:
            subdomain: The SubDomain to transform.
            level: The Level of the transformation.
            domain: The Domain of the transformation.
            **kwargs: Additional keyword arguments for the transformation.

                clip (bool): If True, clips the repeated elements to the domain. Default is False.

        Returns:
            SubDomain: The transformed SubDomain.
        """
        if level in self.excluded_levels:
            return subdomain
        clip = kwargs.get("clip", True)
        clipping_plane = domain.profile if clip else None

        return repeat(
            subdomain,
            self.repeat_info,
            self.repeat_info_2,
            clipping_plane,
            domain.grid_size,
        )


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
        self, subdomain: SubDomain, level: int, domain: Domain, **kwargs: Any
    ) -> SubDomain:
        """Apply the diffeomorphism to the polygon."""
        if isinstance(subdomain, Polygon):
            return self.__map_polygon(subdomain)
        elif isinstance(subdomain, MultiPolygon):
            return MultiPolygon([self.__map_polygon(poly) for poly in subdomain])
        else:
            raise ValueError(f"Unknown polygon type {type(subdomain)}.")

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
