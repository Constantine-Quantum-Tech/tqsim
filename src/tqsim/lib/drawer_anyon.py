# This code is part of TQSim.
#
# (C) Copyright Constantine Quantum Technologies, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from typing import Any
import numpy as np
import numpy.typing as npt


class DrawerAnyon:
    def __init__(self, initial_id: Any, pos: float) -> None:
        self.__i = np.linspace(0, 1)
        self.__one = np.ones_like(self.__i)
        self.__x = [np.copy(self.__i)]
        self.__y = [pos * self.__one]

        self.__color = "black"
        self.__label = ""

    def __repr__(self) -> str:
        return f"{self.__label} (color: {self.__color}"

    def __str__(self) -> str:
        return f"{self.__label} (color: {self.__color}"

    @property
    def label(self) -> str:
        return self.__label

    @label.setter
    def label(self, new: str) -> None:
        self.__label = new

    @property
    def color(self) -> str:
        return self.__color

    @color.setter
    def color(self, new: str) -> None:
        self.__color = new

    @property
    def x(self) -> list[npt.NDArray[np.float64]]:
        return self.__x

    @x.setter
    def x(self, new: npt.NDArray[np.float64]) -> None:
        self.__x.append(new)

    @property
    def y(self) -> list[npt.NDArray[np.float64]]:
        return self.__y

    @y.setter
    def y(self, new: npt.NDArray[np.float64]) -> None:
        self.__y.append(new)

    def get_last_x(self) -> np.floating[Any]:
        return self.__x[-1][-1]

    def get_first_y(self) -> np.floating[Any]:
        return self.__y[0][0]

    def get_last_y(self) -> np.floating[Any]:
        return self.__y[-1][-1]

    def add_identity(self) -> None:
        self.__x.append(self.__i + self.get_last_x())
        self.__y.append(self.__one * self.get_last_y())
