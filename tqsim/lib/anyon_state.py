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

from abc import ABC, abstractmethod
import numpy as np


class AnyonState(ABC):
    """Abstract base class for anyonic states."""

    @abstractmethod
    def is_valid(self, model):
        pass


class StandardAnyonState(AnyonState):
    r"""Anyonic State involving multiple anyons in a linear
    fusion tree (standard basis).

    Example:

        i0  i1  i2  i3
         \ /   /    /
          o1  /    /
           \ /    /
            o2   /
             \  /
              o3

        inputs:
        -------
            inputs   : np.array([i0, i1, i2, ..])
            outcomes : np.array([o1, o2, o3, ..])
    """

    def __init__(self, inputs: np.ndarray, outcomes: np.ndarray):
        self.inputs = inputs
        self.outcomes = outcomes
        assert (
            len(outcomes) == len(inputs) - 1
        ), "Number of outcomes must be one less than number of inputs."

    def __eq__(self, value):
        return np.array_equal(self.inputs, value.inputs) and np.array_equal(
            self.outcomes, value.outcomes
        )

    def is_valid(self, model):
        return model.check_state(self)

    def __repr__(self):
        return (
            f"StandardAnyonState(inputs={self.inputs}, "
            f"outcomes={self.outcomes})"
        )

    def inner_product(self, other):
        """Compute the inner product between two StandardAnyonState instances.

        Parameters
        ----------
        other : StandardAnyonState
            The other state to compute the inner product with.

        Returns
        -------
        float
            1.0 if the states are identical, 0.0 otherwise.
        """
        if not isinstance(other, StandardAnyonState):
            raise ValueError(
                "Inner product (for now) is only defined between states of the same type."
            )

        if np.array_equal(self.inputs, other.inputs) and np.array_equal(
            self.outcomes, other.outcomes
        ):
            return 1.0
        else:
            return 0.0


class SparseAnyonState(AnyonState):
    r"""Anyonic State involving multiple qudits,
    each qudit being made of multiple anyons (sparse basis).

    Example:
        a b ..        g  h
        \/  / \/  / \/  /
        i\ /  k\ /  m\ /
          \     /     /
          j\  l/     /n
            \ /     /
            s\     /
              \   /
               \ /
               t|
        Here,
            charges = [a, b, .., g, h, i, j, k, l, m, n, s, t]
            nb_qudits = 3
            nb_anyons_per_qudit = 3
    """

    def __init__(
        self, charges: np.ndarray, nb_qudits: int, nb_anyons_per_qudit: int
    ):
        self.charges = charges
        self.nb_qudits = nb_qudits
        self.nb_anyons_per_qudit = nb_anyons_per_qudit
        assert len(charges) == nb_qudits * (2 * nb_anyons_per_qudit - 1) + (
            nb_qudits - 1
        ), "Charges length does not match the number of qudits and anyons per qudit."

    def __eq__(self, value):
        return (
            np.array_equal(self.charges, value.charges)
            and self.nb_qudits == value.nb_qudits
            and self.nb_anyons_per_qudit == value.nb_anyons_per_qudit
        )

    def is_valid(self, model):
        return model.check_state(self)

    def __repr__(self):
        # Create a string representation of the SparseAnyonState
        return (
            f"SparseAnyonState(charges={self.charges}, "
            f"nb_qudits={self.nb_qudits}, "
            f"nb_anyons_per_qudit={self.nb_anyons_per_qudit})"
        )

    def __getitem__(self, index):
        return self.charges[index]

    def get_qudit_state(self, qudit_index):
        """Retrieve the charges corresponding to a specific qudit
        as StandardAnyonState.

        Parameters
        ----------
        qudit_index : int: 0 to nb_qudits-1
            The index of the qudit to retrieve.

        Returns
        -------
        StandardAnyonState
            The StandardAnyonState corresponding to the specified qudit.
        """
        assert 0 <= qudit_index < self.nb_qudits, "Invalid qudit index."
        qudit_len = self.nb_anyons_per_qudit - 1
        inputs_start = qudit_index * self.nb_anyons_per_qudit
        inputs_end = inputs_start + self.nb_anyons_per_qudit

        outcomes_start = (self.nb_qudits * self.nb_anyons_per_qudit) + (
            qudit_index * qudit_len
        )
        outcomes_end = outcomes_start + qudit_len

        return StandardAnyonState(
            inputs=self.charges[inputs_start:inputs_end],
            outcomes=self.charges[outcomes_start:outcomes_end],
        )

    def get_inputs(self):
        """Retrieve the inputs of the SparseAnyonState.

        Returns
        -------
        np.ndarray
            The inputs of the SparseAnyonState.
        """
        total_inputs = self.nb_qudits * self.nb_anyons_per_qudit
        return self.charges[:total_inputs]

    def get_outcomes(self):
        """Retrieve the outcomes of the SparseAnyonState.

        Returns
        -------
        np.ndarray
            The outcomes of the SparseAnyonState.
        """
        total_outcomes = self.nb_qudits - 1
        return self.charges[-total_outcomes:]

    def inner_product(self, other):
        """Compute the inner product between two SparseAnyonState instances.

        Parameters
        ----------
        other : SparseAnyonState
            The other state to compute the inner product with.

        Returns
        -------
        float
            1.0 if the states are identical, 0.0 otherwise.
        """
        if not isinstance(other, SparseAnyonState):
            raise ValueError(
                "Inner product (for now) is only defined between states of the same type."
            )

        if (
            np.array_equal(self.charges, other.charges)
            and self.nb_qudits == other.nb_qudits
            and self.nb_anyons_per_qudit == other.nb_anyons_per_qudit
        ):
            return 1.0
        else:
            return 0.0


class ComputationalSparseAnyonState(AnyonState):
    r"""Anyonic State involving multiple qudits,
    each qudit being made of multiple anyons (sparse basis)
    such that the input anyons have identical charges.

    Example:
        a a ..        a  a
        \/  / \/  / \/  /
        i\ /  k\ /  m\ /
          \     /     /
          j\  l/     /n
            \ /     /
            s\     /
              \   /
               \ /
               t|
        Here,
            charges = [i, j, k, l, m, n, s, t]
            nb_qudits = 3
            nb_anyons_per_qudit = 3
    """

    def __init__(
        self,
        charges: np.ndarray,
        nb_qudits: int,
        nb_anyons_per_qudit: int,
        input_charge: int,
    ):
        self.charges = charges
        self.nb_qudits = nb_qudits
        self.nb_anyons_per_qudit = nb_anyons_per_qudit
        self.input_charge = input_charge
        assert len(charges) == nb_qudits * (nb_anyons_per_qudit - 1) + (
            nb_qudits - 1
        ), "Charges length does not match the number of qudits and anyons per qudit."

    def __eq__(self, value):
        return (
            np.array_equal(self.charges, value.charges)
            and self.nb_qudits == value.nb_qudits
            and self.nb_anyons_per_qudit == value.nb_anyons_per_qudit
            and self.input_charge == value.input_charge
        )

    def is_valid(self, model):
        sparse_state = SparseAnyonState(
            charges=np.concatenate(
                (
                    np.full(
                        self.nb_qudits * self.nb_anyons_per_qudit,
                        self.input_charge,
                    ),
                    self.charges,
                )
            ),
            nb_qudits=self.nb_qudits,
            nb_anyons_per_qudit=self.nb_anyons_per_qudit,
        )
        return model.check_state(sparse_state)

    def __repr__(self):
        return (
            f"ComputationalSparseAnyonState(charges={self.charges}, "
            f"nb_qudits={self.nb_qudits}, "
            f"nb_anyons_per_qudit={self.nb_anyons_per_qudit}, "
            f"input_charge={self.input_charge})"
        )

    def inner_product(self, other):
        """Compute the inner product between two ComputationalSparseAnyonState instances.

        Parameters
        ----------
        other : ComputationalSparseAnyonState,SparseAnyonState
            The other state to compute the inner product with.

        Returns
        -------
        float
            1.0 if the states are identical, 0.0 otherwise.
        """
        if isinstance(other, SparseAnyonState):
            sparse_self = SparseAnyonState(
                charges=np.concatenate(
                    (
                        np.full(
                            self.nb_qudits * self.nb_anyons_per_qudit,
                            self.input_charge,
                        ),
                        self.charges,
                    )
                ),
                nb_qudits=self.nb_qudits,
                nb_anyons_per_qudit=self.nb_anyons_per_qudit,
            )
            return sparse_self.inner_product(other)
        elif isinstance(other, ComputationalSparseAnyonState):
            if (
                np.array_equal(self.charges, other.charges)
                and self.nb_qudits == other.nb_qudits
                and self.nb_anyons_per_qudit == other.nb_anyons_per_qudit
                and self.input_charge == other.input_charge
            ):
                return 1.0
            else:
                return 0.0
        else:
            raise ValueError(
                "Inner product (for now) is only defined between states of the same type or with SparseAnyonState."
            )
