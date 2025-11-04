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

import numpy as np

from tqsim.lib.anyon_model import AnyonModel

fusion_matrix = np.zeros((2, 2, 2))
for a1 in range(2):
    for a2 in range(2):
        for outcome in range(2):
            if (a1, a2) == (1, 1):
                fusion_matrix[a1, a2, outcome] = 1
            else:
                if (a1 + a2) == outcome:
                    fusion_matrix[a1, a2, outcome] = 1


def F(a1, a2, a3, outcome):
    """
    F matrix
    """
    inv_phi = (np.sqrt(5) - 1) / 2  # inverse of golden number
    f_matrix = np.array([[0, 0], [0, 0]])

    # a1 + a2 + a3 + outcome = 4
    if a1 + a2 + a3 + outcome == 4:
        f_matrix = np.array([[inv_phi, np.sqrt(inv_phi)], [np.sqrt(inv_phi), -inv_phi]])
    # a1 + a2 + a3 + outcome = 3
    elif a1 + a2 + a3 + outcome == 3:
        f_matrix = np.array([[0, 0], [0, 1]])
    # a1 + a2 + a3 + outcome = 2
    elif a1 + a2 + a3 + outcome == 2:
        if a1 + a2 == 2:
            f_matrix = np.array([[0, 1], [0, 0]])
        elif a2 + a3 == 2:
            f_matrix = np.array([[0, 0], [1, 0]])
        elif a1 + a3 == 2:
            f_matrix = np.array([[0, 0], [0, 1]])
        elif a3 + outcome == 2:
            f_matrix = np.array([[0, 1], [0, 0]])
        elif a1 + outcome == 2:
            f_matrix = np.array([[0, 0], [1, 0]])
        elif a2 + outcome == 2:
            f_matrix = np.array([[0, 0], [0, 1]])
    # a1 + a2 + a3 + outcome = 1
    # a1 + a2 + a3 + outcome = 0
    elif a1 + a2 + a3 + outcome == 0:
        f_matrix = np.array([[1, 0], [0, 0]])

    return f_matrix


def R(a1, a2):
    """
    R matrix
    """
    if a1 + a2 == 2:
        r_matrix = np.array(
            [[np.exp(-4 * np.pi * 1j / 5), 0], [0, np.exp(3 * np.pi * 1j / 5)]]
        )
    else:
        r_matrix = np.array([[1, 0], [0, 1]])

    return r_matrix


F_matrix = np.zeros((2, 2, 2, 2, 2, 2)) * (1 + 0j)
R_matrix = np.zeros((2, 2, 2)) * (1 + 0j)

for a1 in range(2):
    for a2 in range(2):
        for a3 in range(2):
            for outcome in range(2):
                F_matrix[a1, a2, a3, outcome] = F(a1, a2, a3, outcome)

for a1 in range(2):
    for a2 in range(2):
        R_matrix[a1, a2] = R(a1, a2).diagonal()

FIBONACCI_MODEL = AnyonModel(fusion_matrix, F_matrix, R_matrix, name="Fibonacci")
