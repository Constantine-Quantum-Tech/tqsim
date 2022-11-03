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

from typing import List

import numpy as np


def check_rule(anyon1: int, anyon2: int, outcome: int) -> bool:
    """ Returns True if 'anyon1 x anyon2 = outcome' obeys the Fibonacci
    fusion rules, returns False otherwise.

    Parameters
    ----------
    anyon1 : int
        Anyon charge of the 1st anyon.
    anyon2 : int
        Anyon charge of the 2nd anyon.
    outcome : int
        Anyon charge of the fusion result.

    Returns
    -------
    bool
        True if the Fibonacci fusion rules are obeyed, False otherwise.

    """
    if anyon1 and anyon2:
        return True
    elif (anyon1 or anyon2) and outcome == 1:
        return True
    elif not (anyon1 or anyon2) and outcome == 0:
        return True
    else:
        return False


def check_outcomes(outcomes: List[int]) -> bool:
    previous_outcome = 1

    for outcome in outcomes:
        if check_rule(previous_outcome, 1, outcome):
            previous_outcome = outcome
        else:
            return False
    return True


def check_state(state) -> bool:
    nb_qudits = len(state["qudits"])
    qudit_len = len(state["qudits"][0])

    for qudit in state["qudits"]:
        if len(qudit) == qudit_len:
            if not check_outcomes(qudit):
                return False
        else:
            return False

    if nb_qudits != len(state["roots"]) + 1:
        return False

    previous_outcome = state["qudits"][0][-1]

    for i, outcome in enumerate(state["roots"]):
        if check_rule(previous_outcome, state["qudits"][i + 1][-1], outcome):
            previous_outcome = outcome
        else:
            return False
    return True


def gen_state(comb: List[int], nb_qudits: int, qudit_len: int):
    state = {"qudits": [], "roots": []}

    for i, label in enumerate(comb):
        if i < nb_qudits * qudit_len:
            if i % qudit_len:
                state["qudits"][-1].append(label)
            else:
                state["qudits"].append([label])
        else:
            state["roots"].append(label)

    return state


def generate_basis(nb_qudits: int, nb_anyons_per_qudit: int):
    """Generates all the basis states for a system of a given number of
    qudits, and a given number of anyons per qudit.

    Parameters
    ----------
    nb_qudits : int
        Number of qudits in the circuit.
    nb_anyons_per_qudit : int
        Number of anyons in each qudit.

    Returns
    -------
    basis : List[basis states]
        A list of basis states.
    """
    nb_roots = nb_qudits - 1
    qudit_len = nb_anyons_per_qudit - 1
    nb_labels = nb_qudits * qudit_len + nb_roots

    basis = []

    curr_comb = [0] * nb_labels
    final_comb = [1] * nb_labels

    curr_state = gen_state(curr_comb, nb_qudits, qudit_len)

    if check_state(curr_state):
        basis.append(curr_state)

    while not np.all(curr_comb == final_comb):
        for i, label in enumerate(curr_comb):
            if label == 0:
                curr_comb[i] = 1
                break
            else:
                curr_comb[i] = 0

        curr_state = gen_state(curr_comb, nb_qudits, qudit_len)

        if check_state(curr_state):
            basis.append(curr_state)

    return basis
