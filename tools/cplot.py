"""
Created on Thu Aug 27 21:54:22 2020

@author: abduhu

Complex unitary matrix plotting
********


"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from math import pi, sqrt, cos, atan, exp
from colorsys import hls_to_rgb


def cplot(unitary: "np.ndarray", title="", sigma=0.5, show=True, ticks=None):
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
    #sns.set_style("whitegrid", {"grid.linestyle": ":"})

    dims = unitary.shape
    img = []
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
                sat = 0
            else:
                sat = cos(pi * rad / 2) * 0.5 + 0.5
            img[-1].append(hls_to_rgb(hue, lum, sat))

    fig, axs = plt.subplots(nrows=1, figsize=[4, 4])

    if ticks is None:
        plt.xticks([i for i in range(dims[0])])
        plt.yticks([dims[1] - 1 - i for i in range(dims[1])])
    else:
        plt.xticks([i for i in range(dims[0])], labels=ticks)
        plt.yticks(
            [dims[1] - 1 - i for i in range(dims[1])], labels=ticks[::-1]
        )

    plt.grid(linewidth=0.1)
    plt.imshow(img)  # , extent=(1, dims[0], 1, dims[1]))
    axs.xaxis.set_ticks_position("top")
    plt.axis()
    plt.savefig(f"images/{title}.png", dpi=500, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()


def scale(sigma=0.5, title="scale", show=True):
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
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))

    N = 100
    rad = np.tile(np.linspace(0, 1, N).reshape((N, 1)), N)
    theta = np.tile(np.linspace(0, 2*pi, N), (N, 1))

    color = np.ones((N, N, 3))
    #color[:, :, 0] = abs(np.sin(theta/2 + 2*pi/3))
    #color[:, :, 1] = abs(np.sin(theta/2 ))
    #color[:, :, 2] = abs(np.sin(theta/2 + pi/3))
    #color[:, :, 3] = abs(rad)

    for t in range(N):
        for r in range(N):
            hue = t/N
            rdi = sqrt(r/N)
            lum = 0.5 + 0.5 * exp(-rdi / sigma)
            if rdi > 2:
                sat = 0
            else:
                sat = cos(pi * rdi / 2) * 0.5 + 0.5
            color[r, t] = np.array(hls_to_rgb(hue, lum, sat))

    # Plot the color map
    ax.pcolormesh(theta, rad, color)

    # Remove labels and ticks
    ax.set_yticklabels([])
    ax.set_xticklabels(['0', '\u03C0/4', '\u03C0/2', '3\u03C0/4', '\u03C0', '5\u03C0/4', '3\u03C0/2', '7\u03C0/4'])
    plt.grid(linewidth=0.2)
    fig.tight_layout()
    plt.savefig(f"images/{title}.png", dpi=500, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close()