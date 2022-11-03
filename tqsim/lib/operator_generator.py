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

from copy import deepcopy

import numpy as np

from .basis_generator import generate_basis


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


def B(a0, a1, a2, outcome):
    """
    Braiding matrix
    """
    b_matrix = F(a0, a1, a2, outcome) @ R(a1, a2) @ F(a0, a2, a1, outcome).T.conjugate()

    return b_matrix


def sigma(index, state_f, state_i):
    """
    Amplitude of getting state_f by applying the braiding operator
    sigma_{index} on state_i.

    Returns:
        the component (state_f, state_i) of the sigma_{index} matrix
    """
    if index <= 0 or index > len(state_i):
        raise ValueError("index value is not valid!")

    stt_f = [1] + state_f
    stt_i = [1] + state_i

    if index - 2 < 0:
        a0 = 0
    elif index - 2 == 0:
        a0 = 1
    else:
        a0 = state_i[index - 3]

    outcome = state_i[index - 1]
    a = stt_i[index - 1]
    b = stt_f[index - 1]
    amplitude = B(a0, 1, 1, outcome)[a, b]

    ket = stt_i
    ket[index - 1] = b
    bra = stt_f
    if ket == bra:
        braket = 1
    else:
        braket = 0
    return amplitude * braket


def L(k, h, i_, i, jj_, jj):
    """
    L matrix component that is used in calculation of braiding between
    two anyons separated in two qudits.
    (see report)

    Inputs:
        k: int: k
        h: int: i_{m(q-1)}
        i_: int: i'_{mq}
        i: int: i_{mq}
        jj_: list: [i'_{(m+1)1},....i'_{(m+1)q}]
        jj: list: [i_{(m+1)1},....i_{(m+1)q}]
    """
    component = 0 + 0j

    qudit_len = len(jj)
    jjj_ = deepcopy(jj_)
    jjj = deepcopy(jj)
    jjj_ = [1] + jjj_
    jjj = [1] + jjj

    init_p = [0] * qudit_len
    final_p = [1] * qudit_len
    new_p = init_p
    while new_p != final_p:
        pp = deepcopy(new_p)
        pp.append(k)
        product = 1 + 0j
        for ii in range(qudit_len):
            product = (
                product
                * F(i, jjj[ii], 1, pp[ii + 1]).T.conjugate()[jjj[ii + 1], pp[ii]]
                * F(i_, jjj_[ii], 1, pp[ii + 1])[pp[ii], jjj_[ii + 1]]
            )

        product = product * B(h, 1, 1, pp[0])[i, i_]
        component += product
        # iterate
        for ii, label in enumerate(new_p):
            if label == 0:
                new_p[ii] = 1
                break
            else:
                new_p[ii] = 0

    # final iteration
    pp = deepcopy(new_p)
    pp.append(k)
    product = 1 + 0j
    for ii in range(qudit_len):
        product = (
            product
            * F(i, jjj[ii], 1, pp[ii + 1]).T.conjugate()[jjj[ii + 1], pp[ii]]
            * F(i_, jjj_[ii], 1, pp[ii + 1])[pp[ii], jjj_[ii + 1]]
        )

    product = product * B(h, 1, 1, pp[0])[i, i_]
    component += product

    return component


def S(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj):
    """
    S matrix or sewing matrix is used in calculation of braiding operator
    between two anyons separated between two qudits not fused imedialtely.

    Inputs:
        jm: int: j_m
        jmo: int: j_{m-1}
        jmoo: int: j_{m-2}
        jmo_: int: j'_{m-1}
        h: int: i_{m(q-1)}
        i_: int: i'_{mq}
        i: int: i_{mq}
        jj_: list: [i'_{(m+1)1},....i'_{(m+1)q}]
        jj: list: [i_{(m+1)1},....i_{(m+1)q}]
    """
    component = 0 + 0j

    for kk in [0, 1]:
        component += (
            F(jmoo, i, jj[-1], jm)[jmo, kk]
            * L(kk, h, i_, i, jj_, jj)
            * F(jmoo, i_, jj_[-1], jm).T.conjugate()[kk, jmo_]
        )

    return component


def gen_sigma(index, state_i, state_f):
    qudit_len = len(state_i["qudits"][0])
    nb_anyons_per_qudit = qudit_len + 1

    amplitude = 0
    braket = 1

    if index % nb_anyons_per_qudit > 0:

        m = index // nb_anyons_per_qudit
        idx = index % nb_anyons_per_qudit
        amplitude = sigma(idx, state_f["qudits"][m], state_i["qudits"][m])

        for i, qudit in enumerate(state_i["qudits"]):
            if i == m:
                continue
            elif qudit != state_f["qudits"][i]:
                braket = 0

        for i, root in enumerate(state_i["roots"]):
            if root != state_f["roots"][i]:
                braket = 0

    else:
        m = (index // nb_anyons_per_qudit) - 1

        new_state_i = deepcopy(state_i)
        new_state_i["qudits"][m][-1] = deepcopy(state_f["qudits"][m][-1])
        new_state_i["qudits"][m + 1] = deepcopy(state_f["qudits"][m + 1])
        """
            jm: int: j_m
            jmo: int: j_{m-1}
            jmoo: int: j_{m-2}
            jmo_: int: j'_{m-1}
            h: int: i_{m(q-1)}
            i_: int: i'_{mq}
            i: int: i_{mq}
            jj_: list: [i'_{(m+1)1},....i'_{(m+1)q}]
            jj: list: [i_{(m+1)1},....i_{(m+1)q}]
        """
        if m + 1 > 2:
            new_state_i["roots"][m - 1] = state_f["roots"][m - 1]

            jj_ = deepcopy(new_state_i["qudits"][m + 1])
            jj = deepcopy(state_i["qudits"][m + 1])
            h = state_i["qudits"][m][-2]
            i = state_i["qudits"][m][-1]
            i_ = new_state_i["qudits"][m][-1]

            jmo_ = new_state_i["roots"][m - 1]
            jmoo = state_i["roots"][m - 2]
            jmo = state_i["roots"][m - 1]
            jm = state_i["roots"][m]

        elif m + 1 == 2:
            new_state_i["roots"][m - 1] = state_f["roots"][m - 1]

            jj_ = deepcopy(new_state_i["qudits"][m + 1])
            jj = deepcopy(state_i["qudits"][m + 1])
            h = state_i["qudits"][m][-2]
            i = state_i["qudits"][m][-1]
            i_ = new_state_i["qudits"][m][-1]

            jmo_ = new_state_i["roots"][m - 1]
            jmoo = state_i["qudits"][0][-1]
            jmo = state_i["roots"][m - 1]
            jm = state_i["roots"][m]

        elif m + 1 == 1:

            jj_ = deepcopy(new_state_i["qudits"][m + 1])
            jj = deepcopy(state_i["qudits"][m + 1])
            h = state_i["qudits"][m][-2]
            i = state_i["qudits"][m][-1]
            i_ = new_state_i["qudits"][m][-1]

            jmo_ = new_state_i["qudits"][0][-1]
            jmoo = 0
            jmo = state_i["qudits"][0][-1]
            jm = state_i["roots"][m]

        amplitude += S(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        if new_state_i != state_f:
            braket = 0

    return braket * amplitude


def generate_braiding_operator(index: int, nb_qudits: int, nb_anyons_per_qudit: int):
    """Generates the braiding operator of index 'index' for a system of
    a given number of qudits and anyons per qudit.
    This operator braids anyons at positions 'index' and 'index'+1.

    Parameters
    ----------
    index : int
        The operator's index.
    nb_qudits : int
        Number of qudits in the circuit.
    nb_anyons_per_qudit : int
        Number of anyons in each qudit.

    Returns
    -------
    List
        Matrix representation of the braiding operator.

    """
    basis = generate_basis(nb_qudits, nb_anyons_per_qudit)

    sigmas = []
    for f, base_f in enumerate(basis):
        sigmas.append([])
        for base_i in basis:
            sigmas[f].append(gen_sigma(index, base_i, base_f))

    return sigmas
