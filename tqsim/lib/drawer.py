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

from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .drawer_anyon import DrawerAnyon
from .utils import matplotlib_close_if_inline


class Drawer:
    def __init__(
        self,
        nb_qudits: int,
        nb_anyons_per_qudit: int,
    ):
        self._i = np.linspace(0, 1)
        self._one = np.ones_like(self._i)
        self._hi = np.linspace(0, 0.46, 25)
        self._id = 2 * (0.5 - self._hi[-1])

        self.__nb_qudits = nb_qudits
        self.__nb_anyons_per_qudit = nb_anyons_per_qudit
        self.__nb_anyons = nb_qudits * nb_anyons_per_qudit

        self.__anyons = {}
        self.__idx_map = {}

        self.__STARTING_INDEX = 1

        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudit):
                id = i * self.__nb_anyons_per_qudit + j + self.__STARTING_INDEX
                self.__idx_map[id] = id
                self.__anyons[id] = DrawerAnyon(
                    id, (i * (self.__nb_anyons_per_qudit + 1) + j)
                )

    @property
    def anyons(self):
        return self.__anyons

    def _sigmoid(self, x):
        expr = 1 / (1 + np.exp(-(x - 0.5) * 12))
        expr -= expr[0]
        expr /= expr[-1]
        return expr

    def braid(self, m: int, n: int) -> None:

        # Preprocess
        m_init = self.__idx_map[m]
        n_init = self.__idx_map[n]

        anyon_over = self.__anyons[m_init]
        anyon_under = self.__anyons[n_init]
        anyons_idle = map(
            lambda x: self.__anyons[x],
            filter(lambda x: x != m_init and x != n_init, self.__anyons.keys()),
        )

        distance = abs(anyon_over.get_last_y() - anyon_under.get_last_y())

        # Process the idle anyons
        for curr_anyon in anyons_idle:
            curr_anyon.add_identity()

        # Process the over anyon
        anyon_over.x = anyon_over.get_last_x() + self._i
        anyon_over.y = anyon_over.get_last_y() - distance * np.sign(
            m - n
        ) * self._sigmoid(self._i)

        # Process the bottom anyon
        anyon_under.x = anyon_under.get_last_x() + self._hi
        anyon_under.x = anyon_under.get_last_x() + self._id + self._hi

        sigm = self._sigmoid(np.append(self._hi, np.linspace(1 - self._hi[-1], 1, 25)))

        start_y = anyon_under.get_last_y()

        anyon_under.y = start_y + distance * np.sign(m - n) * sigm[:25]

        anyon_under.y = start_y + distance * np.sign(m - n) * sigm[25:]

        # Renaming
        self.__idx_map[n], self.__idx_map[m] = self.__idx_map[m], self.__idx_map[n]

    def __fuse(self, idx_anyon_top, idx_anyon_bot):
        self.__anyons[idx_anyon_bot].x = (
            self.__anyons[idx_anyon_bot].get_last_x() + self._i
        )

        self.__anyons[idx_anyon_top].x = (
            self.__anyons[idx_anyon_top].get_last_x() + self._i
        )

        distance = 0.5 * abs(
            self.__anyons[idx_anyon_bot].get_last_y()
            - self.__anyons[idx_anyon_top].get_last_y()
        )

        self.__anyons[idx_anyon_bot].y = self.__anyons[
            idx_anyon_bot
        ].get_last_y() + distance * self._sigmoid(self._i)

        self.__anyons[idx_anyon_top].y = self.__anyons[
            idx_anyon_top
        ].get_last_y() - distance * self._sigmoid(self._i)

    def measure(self):
        for curr_anyon in self.__anyons.values():
            curr_anyon.add_identity()

        # Fusing anyons within qudits
        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudit - 1):
                # Idle anyons
                for k in range(j + 2, self.__nb_anyons_per_qudit):
                    final_idx = (
                        i * self.__nb_anyons_per_qudit + k + self.__STARTING_INDEX
                    )
                    idx = self.__idx_map[final_idx]
                    self.__anyons[idx].add_identity()

                # Fusing
                final_bot_idx = (
                    i * self.__nb_anyons_per_qudit + j + self.__STARTING_INDEX
                )

                idx_anyon_bot = self.__idx_map[final_bot_idx]
                idx_anyon_top = self.__idx_map[final_bot_idx + 1]

                self.__fuse(idx_anyon_top, idx_anyon_bot)

        # Fusing qudits
        for i in range(1, self.__nb_qudits):
            # 1 -> None
            # 2 -> 1
            # 3 -> 1, 2

            # Idle qudits
            for k in range(i + 1, self.__nb_qudits):
                # 2 -> None
                # 3 -> 2
                final_idx = (
                    (k + 1) * self.__nb_anyons_per_qudit + self.__STARTING_INDEX - 1
                )
                idx = self.__idx_map[final_idx]
                self.__anyons[idx].add_identity()

            # Fusing
            final_bot_idx = i * self.__nb_anyons_per_qudit + self.__STARTING_INDEX - 1
            final_top_idx = (
                (i + 1) * self.__nb_anyons_per_qudit + self.__STARTING_INDEX - 1
            )

            idx_anyon_bot = self.__idx_map[final_bot_idx]
            idx_anyon_top = self.__idx_map[final_top_idx]

            self.__fuse(idx_anyon_top, idx_anyon_bot)

    def draw(self):
        width = self.__anyons[1].get_last_x() * 0.5
        height = self.__nb_anyons * 0.3
        fig, ax = plt.subplots(1, 1, figsize=(width, height))

        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudit):
                k = i * self.__nb_anyons_per_qudit + j + 1
                curr_anyon = self.__anyons[k]

                for x, y in zip(curr_anyon.x, curr_anyon.y):
                    ax.plot(x, y, curr_anyon.color)
                ax.text(
                    -0.2,
                    curr_anyon.get_first_y(),
                    curr_anyon.label,
                    horizontalalignment="right",
                )
        ax.axis("off")
        if fig:
            matplotlib_close_if_inline(fig)
        return fig
