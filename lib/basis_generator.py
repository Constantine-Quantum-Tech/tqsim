"""
Module for generating basis states for anyonic systems.
"""
import numpy as np
from copy import deepcopy
from lib.anyon_state import AnyonState, StandardAnyonState, SparseAnyonState
from lib.anyon_model import AnyonModel



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

        curr_state = StandardAnyonState(curr_comb[:nb_anyons], curr_comb[nb_anyons:])

        if curr_state.is_valid(self.model):
            basis.append(deepcopy(curr_state))

        while not np.all(curr_comb == final_comb):
            # Increment curr_comb as a binary counter using numpy
            idx = np.argmax(curr_comb == 0)
            curr_comb[:idx] = 0
            curr_comb[idx] = 1

            curr_state = StandardAnyonState(curr_comb[:nb_anyons], curr_comb[nb_anyons:])

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

        curr_comb = np.zeros(nb_labels, dtype=int)
        final_comb = np.ones(nb_labels, dtype=int)

        curr_state = SparseAnyonState(curr_comb, nb_qudits, nb_anyons_per_qudit)

        if curr_state.is_valid(self.model):
            basis.append(deepcopy(curr_state))

        while not np.all(curr_comb == final_comb):
            # Increment curr_comb as a binary counter using numpy
            idx = np.argmax(curr_comb == 0)
            curr_comb[:idx] = 0
            curr_comb[idx] = 1

            curr_state = SparseAnyonState(curr_comb, nb_qudits, nb_anyons_per_qudit)

            if curr_state.is_valid(self.model):
                basis.append(deepcopy(curr_state))

        return basis
