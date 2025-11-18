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


def check_rule(anyon_1, anyon_2, outcome):
    r"""
    anyons can be either 0 or 1 or 2
        0 : vacuum
        1 : ising anyon
        2 : ising fermion

    1 x 1 = 0 + 2
    1 X 2 = 1
    2 x 2 = 0
    0 x 0 = 0
    0 x a = a
    """
    check = False
    if anyon_1 == 1 and anyon_2 == 2 or anyon_1 == 2 and anyon_2 == 1:
        if outcome == 1:
            check = True
    elif anyon_1 == 0 or anyon_2 == 0:
        if anyon_1 + anyon_2 == outcome:
            check = True
    elif anyon_1 + anyon_2 == outcome or (anyon_1 + anyon_2) % 2 == outcome:
        check = True

    return check


def f(a1, a2, a3, outcome):
    """
    F matrix for Ising model
    """
    f_matrix = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    if a1 == 1 and a2 == 1 and a3 == 1:
        if outcome == 1:
            return np.array([[1, 0, 1], [0, 0, 0], [1, 0, -1]]) / np.sqrt(2)
        elif outcome == 2:
            return np.array([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
    else:

        possible_i = []
        for ii in [0, 1, 2]:
            if check_rule(a1, a2, ii) and check_rule(ii, a3, outcome):

                possible_i.append(ii)

        possible_j = []
        for jj in [0, 1, 2]:
            if check_rule(a2, a3, jj) and check_rule(a1, jj, outcome):

                possible_j.append(jj)

        if len(possible_i) > 0:
            f_matrix[possible_i[0], possible_j[0]] = 1

    return f_matrix


def r(a1, a2):
    """
    R matrix
    """
    if a1 == 1 and a2 == 1:
        return np.array(
            [
                [np.exp(-np.pi * 1j / 8), 0, 0],
                [0, 0, 0],
                [0, 0, np.exp(3 * np.pi * 1j / 8)],
            ]
        )

    elif a1 == 1 and a2 == 2 or a1 == 2 and a2 == 1:
        return np.array([[0, 0, 0], [0, 1j, 0], [0, 0, 0]])

    elif a1 == 2 and a2 == 2:
        return np.array([[-1, 0, 0], [0, 0, 0], [0, 0, 0]])

    elif a1 == 0 or a2 == 0:
        r_matrix = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

        r_matrix[a1 + a2, a1 + a2] = 1
        return r_matrix

    else:
        return np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])


fusion_matrix = np.zeros((3, 3, 3))
for a1 in range(3):
    for a2 in range(3):
        for outcome in range(3):
            if check_rule(a1, a2, outcome):
                fusion_matrix[a1, a2, outcome] = 1

f_matrix = np.zeros((3, 3, 3, 3, 3, 3)) * (1 + 0j)
r_matrix = np.zeros((3, 3, 3)) * (1 + 0j)

for a1 in range(3):
    for a2 in range(3):
        for a3 in range(3):
            for outcome in range(3):
                f_matrix[a1, a2, a3, outcome] = f(a1, a2, a3, outcome)

for a1 in range(3):
    for a2 in range(3):
        r_matrix[a1, a2] = r(a1, a2).diagonal()

ISING_MODEL = AnyonModel(fusion_matrix, f_matrix, r_matrix, name="Ising")