#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Created on Sun Mar  6 12:34:34 2022

@author: appo

Ising Model
-----------

This script defines Ising anyonic model (SU(2)_2) by their fusion rules and
computes braiding generators of any possible qudit represented by Ising anyons.

This model is designed to:

    + check fusion rule validity.
    + check ayonic state validity.
    + generates the basis of fusion space in the standard form (left to right
                                                                fusion order).
    + defines F and R matrices.
    + generates B (braiding) matrices.
    + calculates braiding generators (sigma_n).
"""

import numpy as np
from copy import deepcopy
from tools.cplot import cplot


def F(a1, a2, a3, outcome):
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


def R(a1, a2):
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


def check_state(outcomes):
    r"""
    checks if a state is valid in Ising models. Ex:

        1 1 1 1
        \/ / /
        i\/ /
         j\/
          k\
        outcomes of the state |((((1, 1)_i, 1)_j, 1)_k, 1)_l| are [i, j, k, l]
    Inputs:
        outcomes: list:
            outcomes of the fusion tree by order (L to R)
    """
    check = True
    previous_outcome = 1
    for outcome in outcomes:
        if check_rule(previous_outcome, 1, outcome):
            previous_outcome = deepcopy(outcome)
        else:
            check = False
            break

    return check


def find_basis(n_anyons):
    r"""
    generates all states that form the basis of Hilbert space of n_anyons.
    Inputs:
        n_anyons: in:
            number of anyons.
    Returns:
        List[List]: list of states with their labeling outcomes.
    """
    n_labels = n_anyons - 1

    # Generate all combinations and check if
    # they verify Fibonacci rules
    # (To do) combinations can be generated with binary methods.

    init_comb = [0] * n_labels
    final_comb = [2] * n_labels
    new_comb = init_comb
    states = []

    if check_state(new_comb):
        new_state = deepcopy(new_comb)
        states.append(new_state)

    while not new_comb == final_comb:
        for i, label in enumerate(new_comb):
            if label in [0, 1]:
                new_comb[i] += 1
                break
            else:
                new_comb[i] = 0

        if check_state(new_comb):
            new_state = deepcopy(new_comb)
            states.append(new_state)

    return states


def iterate(n_labels):
    """ """

    init_comb = [0] * n_labels
    final_comb = [2] * n_labels
    new_comb = init_comb

    yield new_comb
    while not new_comb == final_comb:
        for i, label in enumerate(new_comb):
            if label in [0, 1]:
                new_comb[i] += 1
                break
            else:
                new_comb[i] = 0

        yield new_comb


def B(a0, a1, a2, outcome):
    """
    Braiding matrix
    """
    return F(a0, a1, a2, outcome)\
        @ R(a1, a2) @ F(a0, a2, a1, outcome).conjugate().T


def sigma(index, state_f, state_i):
    r"""
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

    ket = stt_i
    ket[index - 1] = b
    bra = stt_f
    if ket != bra:
        return 0

    return B(a0, 1, 1, outcome)[a, b]


def braiding_generator(index, n_anyons, show=True):
    r"""
    calculates the matrix of the braiding generator that exchange
    index'th anyon with the (index + 1)'th anyon.
    Inputs:
        index: int:
        n_anyons: int:
            number of anyons.
    Returns:
        (numpy.array whose dimension equals to the dimension of
        anyons' Hilbert space, basis)
    """

    # basis of Hilbert space
    basis = find_basis(n_anyons)

    # compute components of the braiding matrix
    sig = []
    for f, state_f in enumerate(basis):
        sig.append([])
        for i, state_i in enumerate(basis):
            sig[f].append(sigma(index, state_f, state_i))
    if show:
        cplot(np.array(sig).astype(complex))

    return sig, basis
