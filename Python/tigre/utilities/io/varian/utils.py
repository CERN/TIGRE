from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
import os
from typing import Callable, Union
import xml.etree.ElementTree as ET

PathLike = Union[str, os.PathLike]
XMLReader = Callable[[PathLike], ET.Element]


class XML:
    def __init__(self, filepath: PathLike, xml_reader: XMLReader) -> None:
        self.name = os.path.basename(filepath)
        self.root = xml_reader(filepath)
        self.ns = self._get_namespace()

    def _get_namespace(self) -> dict[str, str]:
        """Returns the namespace of the xml file."""
        return {"": self.root.tag.split("}")[0].strip("{")}

    def _get_field(self, field: str) -> ET.Element:
        elem = self.root.find(field, self.ns)
        if elem is None:
            basefield = field.split("/")[-1]
            raise RuntimeError(f"{basefield} not found in {self.name}")
        else:
            return elem


def sort_mod_N(a: NDArray, N: int) -> tuple[NDArray, NDArray]:
    a = a % N
    i_sort = np.argsort(a)
    return a, i_sort


def interp_weight(x: NDArray, xp: NDArray, N: int = 360) -> tuple[int, int, NDArray]:
    """Linear interpolation of xp at x, mod(N). xp must be monotonically increasing"""

    i_min = int(np.argmin(abs(x - xp)))
    if x - xp[i_min] >= 0:
        i_lower = i_min
    else:
        i_lower = (i_min - 1) % len(xp)

    i_upper = (i_lower + 1) % len(xp)

    xp_diff = (xp[i_upper] - xp[i_lower]) % N
    weights = (np.array([(xp[i_upper] - x), (x - xp[i_lower])]) % N) / float(xp_diff)
    return i_lower, i_upper, weights


def mm2cm(x: int | float | NDArray) -> float | NDArray:
    return 0.1 * x


def cm2mm(x: int | float | NDArray) -> float | NDArray:
    return 10.0 * x
