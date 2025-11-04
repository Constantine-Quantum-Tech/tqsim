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
# Module for generating basis states for anyonic systems.

import itertools
from copy import deepcopy

import numpy as np

from tqsim.lib.anyon_model import AnyonModel
from tqsim.lib.anyon_state import (
    ComputationalSparseAnyonState,
    SparseAnyonState,
    StandardAnyonState,
)


class BasisGenerator:
    """Abstract class to generate basis states for anyonic systems according to
    the specified model and state structure."""

    pass


class StandardBasisGenerator(BasisGenerator):
    """Generates basis states for a system of anyons in the standard basis."""

    def __init__(self, model: AnyonModel):
        self.model = model

    def generate_basis(self, nb_anyons: int):
        """Generates all the basis states for a system of a given number of anyons.

        Parameters
        ----------
        nb_anyons : int
            The number of anyons in the system.

        Returns
        -------
        List[AnyonState]
            A list of valid AnyonState instances representing the basis states.

        """
        nb_roots = nb_anyons - 1
        nb_labels = nb_anyons + nb_roots

        basis = []

        curr_comb = np.zeros(nb_labels, dtype=int)
        final_comb = np.ones(nb_labels, dtype=int)

        curr_state = StandardAnyonState(
            deepcopy(curr_comb[:nb_anyons]), deepcopy(curr_comb[nb_anyons:])
        )

        if curr_state.is_valid(self.model):
            basis.append(deepcopy(curr_state))

        while not np.all(curr_comb == final_comb):
            # Increment curr_comb as a binary counter using numpy
            idx = np.argmax(curr_comb == 0)
            curr_comb[:idx] = 0
            curr_comb[idx] = 1

            curr_state = StandardAnyonState(
                deepcopy(curr_comb[:nb_anyons]),
                deepcopy(curr_comb[nb_anyons:]),
            )

            if curr_state.is_valid(self.model):
                basis.append(deepcopy(curr_state))

        return basis


class SparseBasisGenerator(BasisGenerator):
    """Generates basis states for a system of anyons in the sparse basis."""

    def __init__(self, model: AnyonModel):
        self.model = model

    def generate_basis(self, nb_qudits: int, nb_anyons_per_qudit: int):
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
        basis : List[SparseAnyonState]
            A list of basis states in the Sparse basis.
        """
        nb_roots = nb_qudits - 1
        qudit_len = nb_anyons_per_qudit - 1
        nb_labels = nb_qudits * (2 * qudit_len + 1) + nb_roots

        basis = []

        for curr_comb in itertools.product(
            [i for i in range(self.model.nb_charges)], repeat=nb_labels
        ):

            curr_state = SparseAnyonState(
                np.array(curr_comb), nb_qudits, nb_anyons_per_qudit
            )

            if curr_state.is_valid(self.model):

                basis.append(deepcopy(curr_state))

        return basis


class ComputationalSparseBasisGenerator(BasisGenerator):
    """Generates Computational basis states for a system of anyons
    in the sparse basis."""

    def __init__(self, model: AnyonModel):
        self.model = model

    def generate_basis(
        self, nb_qudits: int, nb_anyons_per_qudit: int, input_charge: int
    ):
        """Generates all the computational basis states for a system of
        - a given number of qudits,
        - a given number of anyons per qudit, and
        - a specific input anyon charge

        Parameters
        ----------
        nb_qudits : int
            Number of qudits in the circuit.
        nb_anyons_per_qudit : int
            Number of anyons in each qudit.
        input_charge: int
            Charge of the input anyons

        Returns
        -------
        basis : List[SparseAnyonState]
            A list of basis states in the Sparse basis.
        """
        nb_roots = nb_qudits - 1
        qudit_len = nb_anyons_per_qudit - 1
        nb_labels = nb_qudits * (qudit_len) + nb_roots

        basis = []

        for curr_comb in itertools.product(
            [i for i in range(self.model.nb_charges)], repeat=nb_labels
        ):

            curr_state = ComputationalSparseAnyonState(
                np.array(curr_comb),
                nb_qudits,
                nb_anyons_per_qudit,
                input_charge,
            )

            if curr_state.is_valid(self.model):

                basis.append(deepcopy(curr_state))

        return basis
