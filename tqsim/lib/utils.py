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


# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2018.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import matplotlib.pyplot as plt
import matplotlib as mpl
from math import pi, cos, exp, atan


MATPLOTLIB_INLINE_BACKENDS = {
    "module://ipykernel.pylab.backend_inline",
    "module://matplotlib_inline.backend_inline",
    "nbAgg",
}


def matplotlib_close_if_inline(figure):
    """Close the given matplotlib figure if the backend in use draws figures inline.
    If the backend does not draw figures inline, this does nothing.  This function is to prevent
    duplicate images appearing; the inline backends will capture the figure in preparation and
    display it as well, whereas the drawers want to return the figure to be displayed."""
    # This can only called if figure has already been created, so matplotlib must exist.
    import matplotlib.pyplot

    if matplotlib.get_backend() in MATPLOTLIB_INLINE_BACKENDS:
        matplotlib.pyplot.close(figure)


def cplot(cmatrix, sigma=0.5, title=''):
    """Plots a complex-valued matrix with color coding, and a color map.
    'Sigma' controls how much small values are colored. A lower value will
    emphasize small values more.

    Parameters
    ----------
    cmatrix : ndarray
        A complex-valued matrix.
    sigma : float, optional
        Standard deviation squared. The default is 0.5.
    title : str, optional
        Title of the plotted figure. The default is ''.

    Returns
    -------
    None.

    """
    img = []
    for r, row in enumerate(cmatrix):
        img.append([])
        for c in row:
            y = c.imag
            x = c.real

            if x == 0:
                if y > 0:
                    theta = pi/2
                else:
                    theta = -pi/2
            else:
                theta = atan(y/x)
            if x < 0:
                theta += pi
            rad = 1- exp(-(x**2 + y**2)/sigma)
            img[-1].append([cos(theta/2)**2,
                            cos(theta/2 + 2*pi/3)**2,
                            cos(theta/2 - 2*pi/3)**2,
                            rad])

    mpl.rcParams['figure.figsize'] = (10, 10)
    fig, (pl, sc) = plt.subplots(nrows=1, ncols=2, sharex=False)
                                 #figsize=[8, 25])
    sc.imshow(img)
    pl.imshow(scale(sigma=sigma), extent=(-1, 1, -1, 1))
    pl.set_xlabel('Re')
    pl.set_ylabel('Img')
    pl.grid(True)
    plt.title(title)
    plt.show()
    return

def scale(sigma=0.5):
    """
    Plot the scaling spectrum of the complex plane [-1, 1, -i, i]
    """
    img = []
    sc = 50
    for r in range(sc, -sc, -1):
        img.append([])
        for c  in range(-sc, sc, 1):
            y = (r/sc)
            x = (c/sc)
            if x == 0:
                if y > 0:
                    theta = pi/2
                else:
                    theta = -pi/2
            else:
                theta = atan(y/x)
            if x < 0:
                theta += pi
            rad = 1- exp(-(x**2 + y**2)/sigma)
            img[-1].append([cos(theta/2)**2,
                            cos(theta/2 + 2*pi/3)**2,
                            cos(theta/2 - 2*pi/3)**2,
                            rad])

    return img
