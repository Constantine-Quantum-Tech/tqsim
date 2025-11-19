"""

# This code is part of TQSim.
#
# (C) Copyright Constantine Quantum Technologies, 2025.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
#
Created on Thu Aug 27 21:54:22 2020

Complex unitary matrix plotting
********


"""

from colorsys import hls_to_rgb
from math import atan, cos, exp, pi, sqrt

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns  # type: ignore[import-untyped]


def cplot(unitary: "np.ndarray", title: str = "", sigma: float = 0.5, show: bool = True, ticks: list[str] | None = None) -> None:
    """
    Plots complex matrix using chromatic values.
    """
    sns.set_theme(
        context="paper",
        style="whitegrid",
        palette="Spectral",
        font="DejaVu Sans Mono",
        font_scale=1.3,
        color_codes=True,
        rc=None,
    )
    # sns.set_style("whitegrid", {"grid.linestyle": ":"})

    dims = unitary.shape
    img: list[list[tuple[float, float, float]]] = []
    for r, row in enumerate(unitary):
        img.append([])
        for c in row:
            y = c.imag
            x = c.real

            if x == 0:
                if y > 0:
                    theta = pi / 2
                else:
                    theta = -pi / 2
            else:
                theta = atan(y / x)
            if x < 0:
                theta += pi

            hue = theta / (2 * pi)
            rad = sqrt(x**2 + y**2)
            lum = 0.5 + 0.5 * exp(-rad / sigma)
            if rad > 2:
                sat = 0.0
            else:
                sat = cos(pi * rad / 2) * 0.5 + 0.5
            img[-1].append(hls_to_rgb(hue, lum, sat))

    fig, axs = plt.subplots(nrows=1, figsize=[4, 4])

    if ticks is None:
        plt.xticks([i for i in range(dims[0])])
        plt.yticks([dims[1] - 1 - i for i in range(dims[1])])
    else:
        plt.xticks([i for i in range(dims[0])], labels=ticks)
        plt.yticks([dims[1] - 1 - i for i in range(dims[1])], labels=ticks[::-1])

    plt.grid(linewidth=0.1)
    plt.imshow(img)  # , extent=(1, dims[0], 1, dims[1]))
    axs.xaxis.set_ticks_position("top")
    plt.axis()
    plt.savefig(f"images/{title}.png", dpi=500, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def scale(sigma: float = 0.5, title: str = "scale", show: bool = True) -> None:
    """
    Plot the scaling spectrum of the complex plane [-1, 1, -i, i]
    """
    sns.set_theme(
        context="paper",
        style="whitegrid",
        palette="Spectral",
        font="DejaVu Sans Mono",
        font_scale=1.3,
        color_codes=True,
        rc=None,
    )

    # Create figure and axes
    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))

    num_points = 100
    rad = np.tile(np.linspace(0, 1, num_points).reshape((num_points, 1)), num_points)
    theta = np.tile(np.linspace(0, 2 * pi, num_points), (num_points, 1))

    color = np.ones((num_points, num_points, 3), dtype=float)

    for t in range(num_points):
        for r in range(num_points):
            hue = t / num_points
            rdi = sqrt(r / num_points)
            lum = 0.5 + 0.5 * exp(-rdi / sigma)
            if rdi > 2:
                sat = 0.0
            else:
                sat = cos(pi * rdi / 2) * 0.5 + 0.5
            color[r, t] = np.array(hls_to_rgb(hue, lum, sat))

    # Plot the color map
    ax.pcolormesh(theta, rad, color)

    # Remove labels and ticks
    ax.set_yticklabels([])
    ax.set_xticklabels(
        [
            "0",
            "\u03c0/4",
            "\u03c0/2",
            "3\u03c0/4",
            "\u03c0",
            "5\u03c0/4",
            "3\u03c0/2",
            "7\u03c0/4",
        ]
    )
    plt.grid(linewidth=0.2)
    fig.tight_layout()
    plt.savefig(f"images/{title}.png", dpi=500, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()
