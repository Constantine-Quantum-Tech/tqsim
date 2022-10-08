import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
from .drawer_anyon import DrawerAnyon
from .utils import matplotlib_close_if_inline

class Drawer:
    def __init__(
        self,
        nb_qudits: int,
        nb_anyons_per_qudits: int,
    ):
        self._i = np.linspace(0, 1)
        self._one = np.ones_like(self._i)
        self._hi = np.linspace(0, 0.46, 25)
        self._id = 2 * (0.5 - self._hi[-1])

        self.__nb_qudits = nb_qudits
        self.__nb_anyons_per_qudits = nb_anyons_per_qudits
        self.__nb_anyons = nb_qudits * nb_anyons_per_qudits

        self.anyons: dict[int, DrawerAnyon] = {}

        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudits):
                id = i * self.__nb_anyons_per_qudits + j + 1
                self.anyons[id] = DrawerAnyon(
                    id, (i * (self.__nb_anyons_per_qudits + 1) + j)
                )

    def _sigmoid(self, x):
        expr = 1 / (1 + np.exp(-(x - 0.5) * 12))
        expr -= expr[0]
        expr /= expr[-1]
        return expr

    def braid(self, m: int, n: int) -> None:

        # Preprocess
        anyon_over = self.anyons[m]
        anyon_under = self.anyons[n]
        anyons_idle = map(
            lambda x: self.anyons[x],
            filter(lambda x: x != m and x != n, self.anyons.keys()),
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
        self.anyons[n], self.anyons[m] = self.anyons[m], self.anyons[n]

    def measure(self):
        for curr_anyon in self.anyons.values():
            curr_anyon.add_identity()

        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudits - 1):
                # Idle anyons
                for k in range(j + 2, self.__nb_anyons_per_qudits):
                    idx = i * self.__nb_anyons_per_qudits + k + 1
                    self.anyons[idx].add_identity()

                # Fusing
                idx_anyon_bot = i * self.__nb_anyons_per_qudits + j + 1
                idx_anyon_top = idx_anyon_bot + 1

                self.anyons[idx_anyon_bot].x = (
                    self.anyons[idx_anyon_bot].get_last_x() + self._i
                )

                self.anyons[idx_anyon_top].x = (
                    self.anyons[idx_anyon_top].get_last_x() + self._i
                )

                distance = 0.5 * abs(
                    self.anyons[idx_anyon_bot].get_last_y()
                    - self.anyons[idx_anyon_top].get_last_y()
                )

                self.anyons[idx_anyon_bot].y = self.anyons[
                    idx_anyon_bot
                ].get_last_y() + distance * self._sigmoid(self._i)

                self.anyons[idx_anyon_top].y = self.anyons[
                    idx_anyon_top
                ].get_last_y() - distance * self._sigmoid(self._i)

    def draw(self):
        width = self.anyons[1].get_last_x() * 0.5
        height = self.__nb_anyons * 0.3
        fig, ax = plt.subplots(1, 1, figsize=(width, height))

        for i in range(self.__nb_qudits):
            for j in range(self.__nb_anyons_per_qudits):
                k = i * self.__nb_anyons_per_qudits + j + 1
                curr_anyon = self.anyons[k]

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
