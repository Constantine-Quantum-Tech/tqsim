#!/usr/bin/env python3
r"""
Created on Mon Sep 14 02:55:34 2020

@author: Abdellah Tounsi
******

'fib_multi_qudits' module computes elementary braiding generators of Fibonacci
(SU(2)_3) anyonic model.

Elementary braiding generators (sigma_n) are computed for any number of anyons
grouped in qudits. The general form of the state considered in this module is
illustrated by the following example:
    Example:
        1 1 1 1 1 1 1 1 1
        \/  / \/  / \/  /
        i\ /  k\ /  e\ /
          \     /     /
          j\  l/     /f
            \ /     /
            m\     /
              \   /
               \ /
               t|
        =[(|((1, 1)_i, 1)_j| (X) |((1, 1)_k, 1)_l|)_m (X) |((1, 1)_e, 1)_f|]_t
        state is represented by Python dict
        {'qudits': [[i, j], [k, l], [e, f]], 'roots': [m, t]}

TODO:
    - raise ValueError's
    - Translate to Cpp
"""

from copy import deepcopy

import tqsim.tools.braiding_generators.fib_qudit as fibo
from tqsim.tools.braiding_generators.fib_qudit import braiding_matrix, f_matrix
from tqsim.tools.cplot import cplot


def check_state(state):
    r"""
    Verifies if a state of 'n_qudits' qudit of 'qudit_len' number
    of anyons represented in the tonsorial form is acceptable in
    Fibonacci model.

    Inputs:
        state:
            outcomes of fusion in the tonsorial form.
            Example:
                1 1 1 1 1 1
                \/  / \/  /
                i\ /  k\ /
                  \     /
                  j\  l/
                    \ /
                     |
                    m|
                =(|((1, 1)_i, 1)_j| (X) |((1, 1)_k, 1)_l|)_m
                state is represented by Python dict
                {'qudits': [[i, j], [k, l]], 'roots': [m]}
    """
    check = True

    # Check that all qudits are valid
    qudit_len = len(state["qudits"][0])
    for qudit in state["qudits"]:

        if len(qudit) == qudit_len:
            if not fibo.check_state(qudit):
                check = False
        else:
            check = False

    # Check that the outcomes are valid
    n_qudits = len(state["qudits"])
    if n_qudits != len(state["roots"]) + 1:
        check = False

    previous_outcome = state["qudits"][0][-1]
    for ii, outcome in enumerate(state["roots"]):
        if fibo.check_rule(
            previous_outcome, state["qudits"][ii + 1][-1], outcome
        ):
            previous_outcome = outcome
        else:
            check = False
            break

    return check


def _fill_state(labels, n_qudits, qudit_len, n_roots):
    """Helper function to fill a state from a label list."""
    state = {"qudits": [], "roots": []}
    for ii in range(n_qudits):
        state["qudits"].append([])
        for jj in range(qudit_len):
            state["qudits"][-1].append(0)

    ll = 0
    for ii in range(qudit_len):
        for jj in range(n_qudits):
            state["qudits"][jj][ii] = labels[ll]
            ll += 1

    for ii in range(n_roots):
        state["roots"].append(labels[ll])
        ll += 1

    return state


def _initialize_combinations(n_labels):
    """Initialize starting and ending combinations."""
    new_comb = [0] * n_labels
    final_comb = [1] * n_labels
    return new_comb, final_comb


def _increment_combination(comb):
    """Increment a binary combination in place."""
    for ii, label in enumerate(comb):
        if label == 0:
            comb[ii] = 1
            break
        else:
            comb[ii] = 0


def _generate_valid_states(new_comb, final_comb, n_qudits, qudit_len, n_roots):
    """Generate all valid states from combinations."""
    states = []
    new_state = _fill_state(new_comb, n_qudits, qudit_len, n_roots)
    if check_state(new_state):
        states.append(new_state)

    while new_comb != final_comb:
        _increment_combination(new_comb)
        new_state = _fill_state(new_comb, n_qudits, qudit_len, n_roots)
        if check_state(new_state):
            states.append(new_state)

    return states


def find_basis(n_qudits, qudit_len):
    """
    generates all states that form the basis of Hilbert space
    of anyons grouped by qudits and fused qudit by qudit.

    Inputs:
        n_qudits: int:
            number of qudits.
        qudit_len: int:
            number of outcomes representing one qudit.
    """
    n_roots = n_qudits - 1
    n_labels = n_qudits * qudit_len + n_roots

    new_comb, final_comb = _initialize_combinations(n_labels)
    return _generate_valid_states(
        new_comb, final_comb, n_qudits, qudit_len, n_roots
    )


def find_basis_(n_qudits, qudit_len):
    """
    generates all states that form the basis of Hilbert space
    of anyons grouped by qudits and fused qudit by qudit.

    Inputs:
        n_qudits: int:
            number of qudits.
        qudit_len: int:
            number of outcomes representing one qudit.
    """

    n_roots = n_qudits - 1
    # n_labels = n_qudits * qudit_len + n_roots
    n_anyons_per_qudit = qudit_len + 1

    # generate all combinations and verify if it is valid state
    one_qudit_basis = fibo.find_basis(n_anyons_per_qudit)
    qudit_basis_len = len(one_qudit_basis)

    # iterate roots
    new_comb_roots = []
    final_comb_roots = []
    for _ in range(n_roots):
        new_comb_roots.append(0)
        final_comb_roots.append(1)

    states = []

    # iterate qudits
    new_comb_qudits = []
    final_comb_qudits = []
    for _ in range(n_qudits):
        new_comb_qudits.append(0)
        final_comb_qudits.append(qudit_basis_len - 1)

    qudits = []
    for i in new_comb_qudits:
        qudits.append(one_qudit_basis[i])

    new_state = {}
    new_state["qudits"] = deepcopy(qudits)
    new_state["roots"] = deepcopy(new_comb_roots)

    if check_state(new_state):
        states.append(new_state)

    while not new_comb_qudits == final_comb_qudits:

        for ii, label in enumerate(new_comb_qudits):
            if label < qudit_basis_len - 1:
                new_comb_qudits[ii] += 1
                break
            else:
                new_comb_qudits[ii] = 0

        print(new_comb_qudits)
        qudits = []
        for i in new_comb_qudits:
            qudits.append(one_qudit_basis[i])

        new_state = {}
        new_state["qudits"] = deepcopy(qudits)
        new_state["roots"] = deepcopy(new_comb_roots)
        # print(new_state)

        if check_state(new_state):
            states.append(new_state)

    while not new_comb_roots == final_comb_roots:

        for ii, label in enumerate(new_comb_roots):
            if label == 0:
                new_comb_roots[ii] = 1
                break
            else:
                new_comb_roots[ii] = 0

        # iterate qudits
        new_comb_qudits = []
        final_comb_qudits = []
        for _ in range(n_qudits):
            new_comb_qudits.append(0)
            final_comb_qudits.append(qudit_basis_len - 1)

        qudits = []
        for i in new_comb_qudits:
            qudits.append(one_qudit_basis[i])

        new_state = {}
        new_state["qudits"] = deepcopy(qudits)
        new_state["roots"] = deepcopy(new_comb_roots)
        # print(new_state)

        if check_state(new_state):
            states.append(new_state)

        while not new_comb_qudits == final_comb_qudits:

            for ii, label in enumerate(new_comb_qudits):
                if label < qudit_basis_len - 1:
                    new_comb_qudits[ii] += 1
                    break
                else:
                    new_comb_qudits[ii] = 0

            qudits = []
            for i in new_comb_qudits:
                qudits.append(one_qudit_basis[i])

            new_state = {}
            new_state["qudits"] = deepcopy(qudits)
            new_state["roots"] = deepcopy(new_comb_roots)
            # print(new_state)

            if check_state(new_state):
                states.append(new_state)

    return states


def l_matrix(k, h, i_, i, jj_, jj):
    r"""
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
                * f_matrix(i, jjj[ii], 1, pp[ii + 1]).T.conjugate()[
                    jjj[ii + 1], pp[ii]
                ]
                * f_matrix(i_, jjj_[ii], 1, pp[ii + 1])[pp[ii], jjj_[ii + 1]]
            )

        product = product * braiding_matrix(h, 1, 1, pp[0])[i, i_]
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
            * f_matrix(i, jjj[ii], 1, pp[ii + 1]).T.conjugate()[
                jjj[ii + 1], pp[ii]
            ]
            * f_matrix(i_, jjj_[ii], 1, pp[ii + 1])[pp[ii], jjj_[ii + 1]]
        )

    product = product * braiding_matrix(h, 1, 1, pp[0])[i, i_]
    component += product

    return component


def knitting_matrix(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj):
    r"""
    K matrix or sewing matrix is used in calculation of braiding operator
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
            f_matrix(jmoo, i, jj[-1], jm)[jmo, kk]
            * l_matrix(kk, h, i_, i, jj_, jj)
            * f_matrix(jmoo, i_, jj_[-1], jm).T.conjugate()[kk, jmo_]
        )

    return component


def _validate_sigma_states(state_f_, state_i_):
    """Validate states for sigma computation."""
    if not (check_state(state_f_) or check_state(state_i_)):
        raise ValueError("States are not valid!")


def _check_unchanged_qudits(state_i_, state_f_, m):
    """Check if all qudits except m are unchanged."""
    for ii, qudit in enumerate(state_i_["qudits"]):
        if ii == m:
            continue
        elif qudit != state_f_["qudits"][ii]:
            return False
    return True


def _check_unchanged_roots(state_i_, state_f_):
    """Check if all roots are unchanged."""
    for ii, root in enumerate(state_i_["roots"]):
        if root != state_f_["roots"][ii]:
            return False
    return True


def _compute_within_qudit_sigma(index_, state_f_, state_i_, m):
    """Compute sigma amplitude for braiding within a qudit."""
    amplitude = fibo.sigma(
        index=index_,
        state_f=state_f_["qudits"][m],
        state_i=state_i_["qudits"][m],
    )

    if not _check_unchanged_qudits(state_i_, state_f_, m):
        return 0

    if not _check_unchanged_roots(state_i_, state_f_):
        return 0

    return amplitude


def _prepare_new_state(state_i_, state_f_, m):
    """Prepare the new state for between-qudit braiding."""
    new_state_i = deepcopy(state_i_)
    new_state_i["qudits"][m][-1] = deepcopy(state_f_["qudits"][m][-1])
    new_state_i["qudits"][m + 1] = deepcopy(state_f_["qudits"][m + 1])
    return new_state_i


def _extract_knitting_params_case1(new_state_i, state_i_, state_f_, m):
    """Extract parameters for knitting matrix when m + 1 > 2."""
    new_state_i["roots"][m - 1] = state_f_["roots"][m - 1]
    if new_state_i != state_f_:
        return None

    jj_ = deepcopy(new_state_i["qudits"][m + 1])
    jj = deepcopy(state_i_["qudits"][m + 1])
    h = state_i_["qudits"][m][-2]
    i = state_i_["qudits"][m][-1]
    i_ = new_state_i["qudits"][m][-1]

    jmo_ = new_state_i["roots"][m - 1]
    jmoo = state_i_["roots"][m - 2]
    jmo = state_i_["roots"][m - 1]
    jm = state_i_["roots"][m]

    return (jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)


def _extract_knitting_params_case2(new_state_i, state_i_, state_f_, m):
    """Extract parameters for knitting matrix when m + 1 == 2."""
    new_state_i["roots"][m - 1] = state_f_["roots"][m - 1]
    if new_state_i != state_f_:
        return None

    jj_ = deepcopy(new_state_i["qudits"][m + 1])
    jj = deepcopy(state_i_["qudits"][m + 1])
    h = state_i_["qudits"][m][-2]
    i = state_i_["qudits"][m][-1]
    i_ = new_state_i["qudits"][m][-1]

    jmo_ = new_state_i["roots"][m - 1]
    jmoo = state_i_["qudits"][0][-1]
    jmo = state_i_["roots"][m - 1]
    jm = state_i_["roots"][m]

    return (jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)


def _extract_knitting_params_case3(new_state_i, state_i_, state_f_, m):
    """Extract parameters for knitting matrix when m + 1 == 1."""
    if new_state_i != state_f_:
        return None

    jj_ = deepcopy(new_state_i["qudits"][m + 1])
    jj = deepcopy(state_i_["qudits"][m + 1])
    h = state_i_["qudits"][m][-2]
    i = state_i_["qudits"][m][-1]
    i_ = new_state_i["qudits"][m][-1]

    jmo_ = new_state_i["qudits"][0][-1]
    jmoo = 0
    jmo = state_i_["qudits"][0][-1]
    jm = state_i_["roots"][m]

    return (jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)


def _compute_between_qudits_sigma(state_i_, state_f_, m):
    """Compute sigma amplitude for braiding between qudits."""
    new_state_i = _prepare_new_state(state_i_, state_f_, m)

    if m + 1 > 2:
        params = _extract_knitting_params_case1(new_state_i, state_i_, state_f_, m)
    elif m + 1 == 2:
        params = _extract_knitting_params_case2(new_state_i, state_i_, state_f_, m)
    else:  # m + 1 == 1
        params = _extract_knitting_params_case3(new_state_i, state_i_, state_f_, m)

    if params is None:
        return 0

    return knitting_matrix(*params)


def sigma(index_, state_f_, state_i_):
    """
    Amplitude of getting state_f by applying the braiding operator
    sigma_{index} on state_i.

    Returns:
        the component (state_f, state_i) of the sigma_{index} matrix
    """
    _validate_sigma_states(state_f_, state_i_)

    qudit_len = len(state_i_["qudits"][0])

    # n modulo q > 0: braiding within a qudit
    if index_ % (qudit_len + 1) > 0:
        m = index_ // (qudit_len + 1)
        return _compute_within_qudit_sigma(
            index_ % (qudit_len + 1), state_f_, state_i_, m
        )

    # n modulo q = 0: braiding between qudits
    m = (index_ // (qudit_len + 1)) - 1
    return _compute_between_qudits_sigma(state_i_, state_f_, m)


def braiding_generator(index, n_qudits, qudit_len, show=True):
    """
    calculates matrix representation of the braiding generator -in the basis
    of multi-qudit fusion space- which exchanges
    index'th anyon with the (index + 1)'th anyon.

    Inputs:
        index: int:
            index of braiding operator.
        n_qudits: int:
            number of qudits.
        qudit_len: int:
            number of outcomes representing one qudit.
    Returns:
        (numpy.array whose dimension equals to the dimension of
        anyons' Hilbert space, basis)
    """

    # basis of Hilbert space
    basis = find_basis(n_qudits, qudit_len)

    # compute components of the braiding matrix
    sig = []
    for f, state_f in enumerate(basis):
        sig.append([])
        for i, state_i in enumerate(basis):
            sig[f].append(sigma(index, state_f, state_i))
    if show:
        cplot(sig)

    return sig, basis
