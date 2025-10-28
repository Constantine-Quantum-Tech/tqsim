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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import pi, cos, exp, atan
from typing import List, Tuple, Dict


def einsum_with_names(
    terms: List[Tuple[np.ndarray, Tuple[str, ...]]],
    output_labels: Tuple[str, ...],
) -> np.ndarray:
    """
    terms: list of (array, labels_tuple) where labels_tuple are composite label strings
           e.g. ("i(m,q)","i(m+1,0)","a(m+1,1)","p(3)")
    output_labels: tuple of composite labels for the result (order matters)
    Returns: np.einsum result using the integer-index API so label names can be arbitrary strings.
    """
    # collect all unique labels and assign integer ids
    label_to_int: Dict[str, int] = {}
    next_int = 0

    def get_int(lbl):
        nonlocal next_int
        if lbl not in label_to_int:
            label_to_int[lbl] = next_int
            next_int += 1
        return label_to_int[lbl]

    arrays = []
    index_lists = []

    # convert each term's label tuple into a list of integers
    for arr, labels in terms:
        idxs = [get_int(lbl) for lbl in labels]
        arrays.append(arr)
        index_lists.append(idxs)

    # ensure output labels are in the mapping
    out_idxs = [get_int(lbl) for lbl in output_labels]

    # build the einsum call: np.einsum(arr0, idx0, arr1, idx1, ..., out_idx_list)
    einsum_args = []
    for arr, idxs in zip(arrays, index_lists):
        einsum_args.append(arr)
        einsum_args.append(idxs)
    einsum_args.append(out_idxs)

    return np.einsum(*einsum_args)


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
