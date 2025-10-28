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
from typing import Tuple, List
from copy import deepcopy
from tqsim.lib.utils import einsum_with_names
from tqsim.lib.anyon_state import (
    AnyonState, StandardAnyonState, SparseAnyonState,
    ComputationalSparseAnyonState
)


class AnyonModel:

    def __init__(self, N_symbols, F_matrix, R_matrix, name=None):
        self.N_symbols = N_symbols
        self.F_matrix = F_matrix
        self.R_matrix = R_matrix
        self.braiding_matrix = self._compute_braiding_matrix()

        if name is None:
            self.name == f"model-{np.random.randint(1000)}"
        else:
            self.name = name

    def check_rule(
        self, anyon1: np.ndarray, anyon2: np.ndarray, outcome: np.ndarray
    ) -> np.ndarray:
        """Returns True if 'anyon1 x anyon2 = outcome' obeys the Fibonacci
        fusion rules, returns False otherwise.

        Parameters
        ----------
        anyon1 : array of charges
            Anyon charge of the 1st anyon.
        anyon2 : array of charges
            Anyon charge of the 2nd anyon.
        outcome : array of charges
            Anyon charge of the fusion result.

        Returns
        -------
        array of bool
            True if the Fibonacci fusion rules are obeyed, False otherwise.

        """
        # check that N_symbols[anyon1[i], anyon2[i], outcome[i]] == 1 for all i
        return self.N_symbols[anyon1, anyon2, outcome] == np.ones_like(anyon1)

    def _compute_braiding_matrix(self):
        r"""Computes the braiding matrix for the anyon model.

        [ B_{abc}^j ]_{im} = sum_l [ F_{abc}^j ]_{il} R_{bc}^l
                             [ F_{acb}^j ]^dag_{lm}

        Returns
        -------
        np.ndarray
            The braiding matrix.

        """
        b_matrix = np.einsum(
            "abcjil, bcl, acbjml -> abcjim",
            self.F_matrix,
            self.R_matrix,
            self.F_matrix.conjugate(),
        )
        return b_matrix

    def _compute_L_matrix(self, q: int):
        r"""
        [L_{ a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{p(q+1)}]^{
        i(m,q) i(m+1,1) ... i(m+1,q)
        }_{
        i'(m,q), i'(m+1,1) ... i'(m+1,q)
        }
        =
        \sum_{p(1), .., p(q)}
        prod_{r=1}^{q} [ F_{i(m, q), i(m+1, q-r),
                         a(m+1, q-r+1)}^{p(q-r+2)}
                        ].dag^{i(m+1, q-r+1)}_{p(q-r+1)}
        [ B^{p(1)}_{i(m,q-1), a(m, q), a(m+1, 0)} ]^{i(m, q)}_{i'(m, q)}
        prod_{r=1}^{q} [ F_{i'(m, q), i'(m+1, q-r),
                         a(m+1, q-r+1)}^{p(q-r+2)}
                        ]^{p(q-r+1)}_{i'(m+1, q-r+1)}

        such that
            p_{q+1} = k,
            i'(m+1,0) = i(m+1,0) = a_{m, q}.

        Inputs
        ------
        q : int
            Number of anyons per qudit minus one.
        Returns
        -------
        np.ndarray
            The L matrix.
        """
        assert q > 0, (
            "q must be strictly positive. "
            "For q=1, L is just the braiding matrix."
        )
        terms = []

        # --- Left product of q dagger-F factors
        r = 1
        while r <= q:
            labels = (
                f"i(m,{q})",
                f"i(m+1,{q-r})",
                f"a(m+1,{q-r+1})",
                f"p({q-r+2})",
                f"i(m+1,{q-r+1})",
                f"p({q-r+1})",
            )
            if r == q:
                labels = (
                    f"i(m,{q})",
                    f"a(m,{q})",
                    f"a(m+1,{q-r+1})",
                    f"p({q-r+2})",
                    f"i(m+1,{q-r+1})",
                    f"p({q-r+1})",
                )
            # conj + swap last two axes for dagger
            F_dag = np.conjugate(self.F_matrix).swapaxes(-1, -2)
            terms.append((F_dag, labels))
            r += 1

        # --- B tensor
        B_labels = (
            f"i(m,{q-1})",
            f"a(m,{q})",
            "a(m+1,0)",
            "p(1)",
            f"i(m,{q})",
            f"ip(m,{q})",
        )
        terms.append((self.braiding_matrix, B_labels))

        # --- Right product of q F factors
        r = 1
        while r <= q:
            labels = (
                f"ip(m,{q})",
                f"ip(m+1,{q-r})",
                f"a(m+1,{q-r+1})",
                f"p({q-r+2})",
                f"p({q-r+1})",
                f"ip(m+1,{q-r+1})",
            )
            if r == q:
                labels = (
                    f"ip(m,{q})",
                    f"a(m,{q})",
                    f"a(m+1,{q-r+1})",
                    f"p({q-r+2})",
                    f"p({q-r+1})",
                    f"ip(m+1,{q-r+1})",
                )
            terms.append((self.F_matrix, labels))
            r += 1

        # --- Output indices
        # [L_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{p(q+1)}]^{
        # i(m,q) i(m+1,0) ... i(m+1,q)
        # }_{
        # i'(m,q), i'(m+1,0) ... i'(m+1,q)
        # }
        out_labels = (
            [f"a(m,{q})"]
            + ["a(m+1,0)"]
            + [f"a(m+1,{q-r+1})" for r in range(1, q + 1)]
            + [f"i(m,{q-1})"]
            + [f"p({q+1})"]
            + [f"i(m,{q})"]
            + [f"i(m+1,{r})" for r in range(1, q + 1)]
            + [f"ip(m,{q})"]
            + [f"ip(m+1,{r})" for r in range(1, q + 1)]
        )

        # Call the helper from earlier
        return einsum_with_names(terms, out_labels)

    def compute_knitting_matrix(self, q: int, return_L=False):
        r"""
        See Appendix of https://arxiv.org/abs/2307.01892

        [K_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{j(m-2), j(m)}]^{
        j(m-1) i(m,q)i(m+1,1) ... i(m+1,q)
        }_{
        j'(m-1), i'(m,q), i'(m+1,1) ... i'(m+1,q)
        }
        =
        \sum_{k}
        [ F^{j(m)}_{j(m-2), i(m,q), i(m+1,q)} ]^{j(m-1)}_{k}
        [L_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{k}]^{
        i(m,q) i(m+1,1) ... i(m+1,q)
        }_{
        i'(m,q), i'(m+1,1) ... i'(m+1,q)
        }
        [ F^{j(m)}_{j(m-2), i'(m,q), i'(m+1,q)} ].dagger^{k}_{j'(m-1)}

        Inputs
        ------
        q : int
            Number of anyons per qudit minus one.

        Returns
        -------
        np.ndarray
            The knitting matrix.
            Output labels :

            f"a(m,{q})",
            *[f"a(m+1,{r})" for r in range(0, q + 1)],
            f"i(m,{q-1})",
            "j(m-2)",
            "j(m)",
            "j(m-1)",
            f"i(m,{q})",
            *[f"i(m+1,{r})" for r in range(1, q + 1)],
            "jp(m-1)",
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(1, q + 1)],

        """
        assert q > 0, (
            "q must be strictly positive. "
            "For q=0, K is just the braiding matrix."
        )
        terms = []

        # --- Left F factor ---
        # Example label order (must match how F is actually stored):
        F_labels = (
            "j(m-2)",
            f"i(m,{q})",
            f"i(m+1,{q})",
            "j(m)",
            "j(m-1)",
            "k",
        )
        terms.append((self.F_matrix, F_labels))

        # --- L factor ---
        # L already carries many indices (a, i, i', ...).
        # Here we assume L_labels is known / fixed.
        # Example layout (you must adapt to your actual storage order!):
        # out_labels = (
        #     [f"a(m,{q})"] + [f"a(m+1,{r})" for r in range(0, q+1)] +
        #     [f"i(m,{q-1})"] + [f"p({q+1})"] +
        #     [f"i(m,{q})"] + [f"i(m+1,{r})" for r in range(0, q+1)] +
        #     [f"ip(m,{q})"] + [f"ip(m+1,{r})" for r in range(0, q+1)]
        # )
        L_labels = (
            f"a(m,{q})",
            *[f"a(m+1,{r})" for r in range(0, q + 1)],
            f"i(m,{q-1})",
            "k",
            f"i(m,{q})",
            *[f"i(m+1,{r})" for r in range(1, q + 1)],
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(1, q + 1)],
        )
        L_matrix = self._compute_L_matrix(q)
        terms.append((L_matrix, L_labels))

        # --- Right F dagger factor ---
        # dagger = conjugate and swapaxes(5,4): swap 'k' and 'j(m-1)'
        F_dag = np.conjugate(np.swapaxes(self.F_matrix, 5, 4))
        F_labels_dag = (
            "j(m-2)",
            f"ip(m,{q})",
            f"ip(m+1,{q})",
            "j(m)",
            "k",
            "jp(m-1)",
        )
        terms.append((F_dag, F_labels_dag))

        # --- Output labels for K ---
        # [K_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{j(m-2), j(m)}]^{
        # j(m-1) i(m,q)i(m+1,0) ... i(m+1,q)
        # }_{
        # j'(m-1), i'(m,q), i'(m+1,0) ... i'(m+1,q)
        # }
        out_labels = (
            f"a(m,{q})",
            *[f"a(m+1,{r})" for r in range(0, q + 1)],
            f"i(m,{q-1})",
            "j(m-2)",
            "j(m)",
            "j(m-1)",
            f"i(m,{q})",
            *[f"i(m+1,{r})" for r in range(1, q + 1)],
            "jp(m-1)",
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(1, q + 1)],
        )

        # Perform contraction
        K = einsum_with_names(terms, out_labels)
        if return_L:
            return K, L_matrix
        return K

    def check_state(self, state) -> bool:
        if isinstance(state, SparseAnyonState):
            return self.check_sparse_state(
                state.charges, state.nb_qudits, state.nb_anyons_per_qudit
            )
        elif isinstance(state, StandardAnyonState):
            return self.check_standard_basis_state(
                state.inputs, state.outcomes
            )
        else:
            raise ValueError(
                "State must be either SparseAnyonState or StandardAnyonState"
            )

    def check_standard_basis_state(
        self, inputs: np.ndarray, outcomes: np.ndarray
    ) -> bool:
        r"""
        Check fusion rules validity of a state in the standard basis:

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

        return
        ------
            np.all(self.check_rule(anyon_1, anyon_2, outcomes))

        """
        assert (
            inputs.shape[0] == outcomes.shape[0] + 1
        ), "inputs must have one more element than outcomes"
        anyon_1 = inputs[1::]
        anyon_2 = np.zeros(len(anyon_1), dtype=int)
        anyon_2[0] = inputs[0]
        anyon_2[1::] = outcomes[0:-1]
        rules = self.check_rule(anyon_1, anyon_2, outcomes)
        return np.all(rules)

    def check_sparse_state(
        self, charges: np.ndarray, nb_qudits: int, nb_anyons_per_qudit: int
    ) -> bool:
        r"""
        Check fusion rules validity of a state in the standard basis:

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

        inputs:
        -------
            charges   : np.array of charges per layer as shown above.
            nb_qudits : number of qudits
            nb_anyons_per_qudit : number of anyons per qudit

        return
        ------
            bool

        """
        for i in range(nb_qudits):
            start = i * nb_anyons_per_qudit
            end = (i + 1) * nb_anyons_per_qudit
            inputs = charges[start:end]
            start = nb_qudits * nb_anyons_per_qudit + i * (
                nb_anyons_per_qudit - 1
            )
            end = nb_qudits * nb_anyons_per_qudit + (i + 1) * (
                nb_anyons_per_qudit - 1
            )
            outcomes = charges[start:end]
            check = self.check_standard_basis_state(inputs, outcomes)
            if not check:
                return False

        if nb_qudits == 1:
            return True

        indices = [
            nb_qudits * nb_anyons_per_qudit
            + (i + 1) * (nb_anyons_per_qudit - 1)
            - 1
            for i in range(nb_qudits)
        ]
        inputs = charges[indices]
        outcomes = charges[nb_qudits * (2 * nb_anyons_per_qudit - 1) : :]
        return self.check_standard_basis_state(inputs, outcomes)

    def check_general_state(self) -> bool:
        r"""
        Check fusion rules validity of a state in the general basis:
        Example:
        a b ..        g  h
        \ \/ / \ / /  \ /
         \i /   j /    k
          l/     m    /
           n      \  /
            \      \/
             \     p
              \   /
               \ /
                t
        Here,
            charges = [a, b, .., g, h, i, j, k, l, m, n, t]
            fusion_tree = [
                (g, h, k), .... (m, k, p), (b, .., i), (a, i, l), (l, .., n), (n, p, t)
            ]
            fusion_tree = []
        """
        pass

    def compute_standard_braid_component(
        self,
        initial_state: StandardAnyonState,
        braid_index: int,
        final_state: StandardAnyonState,
    ):
        """
        Computes the probability amplitudes of getting a final anyon state
        from braiding two anyons (of indices i and i+1) in the initial anyon state.

        Inputs:
        -------
            initial_state: StandardAnyonState
            braid_index  : int = i+1 (the index of the first braided anyon)
            final_state  : StandardAnyonState

        Returns:
            Complex number.
        """
        assert braid_index > 0, "braid_index must be greater than 0"
        assert braid_index < len(
            initial_state.inputs
        ), "braid_index must be less than the total number of anyons"
        assert isinstance(
            initial_state, StandardAnyonState
        ), "initial_state must be a StandardAnyonState"
        assert initial_state.is_valid(self), "initial_state must be valid"
        assert final_state.is_valid(self), "final_state must be valid"
        assert len(initial_state.inputs) == len(
            final_state.inputs
        ), "initial_state and final_state must have the same number of anyons"

        intial_inputs = deepcopy(initial_state.inputs)
        final_inputs = deepcopy(final_state.inputs)

        # Permute the anyons i and i+1 in the initial inputs
        temp = intial_inputs[braid_index]
        intial_inputs[braid_index] = intial_inputs[braid_index - 1]
        intial_inputs[braid_index - 1] = temp
        # Check if the permuted initial inputs match the final inputs
        if not np.array_equal(intial_inputs, final_inputs):
            return 0.0 + 0.0j

        initial_outcomes = initial_state.outcomes
        final_outcomes = deepcopy(final_state.outcomes)
        final_outcomes[braid_index - 1] = initial_outcomes[braid_index - 1]
        if not np.array_equal(initial_outcomes, final_outcomes):
            return 0.0 + 0.0j

        a = initial_state.inputs[braid_index - 1]
        b = initial_state.inputs[braid_index]
        if braid_index == len(initial_state.inputs) - 1:
            c = 0  # vacuum
        else:
            c = initial_state.outcomes[braid_index - 1]
        i = initial_state.outcomes[braid_index - 1]
        if braid_index == 1:
            j = 0  # vacuum
        else:
            j = initial_state.outcomes[braid_index - 2]
        m = final_state.outcomes[braid_index - 1]
        amplitude = self.braiding_matrix[a, b, c, j, i, m]
        return amplitude

    def compute_sparse_braid_inner_product(
        self,
        initial_state: SparseAnyonState,
        braid_index: int,
        final_state: SparseAnyonState,
    ):
        """
        Computes the probability amplitudes of getting a final anyon state
        from braiding two anyons (of indices i and i+1) in the initial anyon state.

        Inputs:
        -------
            initial_state: SparseAnyonState
            braid_index  : int = i+1 (the index of the first braided anyon)
            final_state  : SparseAnyonState

        Returns:
            Complex number.
        """
        assert braid_index > 0, "braid_index must be greater than 0"
        assert braid_index < (
            initial_state.nb_qudits * initial_state.nb_anyons_per_qudit
        ), "braid_index must be less than the total number of anyons"
        assert isinstance(
            initial_state, SparseAnyonState
        ), "initial_state must be a SparseAnyonState"
        assert initial_state.is_valid(self), "initial_state must be valid"
        assert final_state.is_valid(self), "final_state must be valid"
        assert (
            initial_state.nb_qudits == final_state.nb_qudits
        ), "initial_state and final_state must have the same number of qudits"
        assert (
            initial_state.nb_anyons_per_qudit
            == final_state.nb_anyons_per_qudit
        ), "initial_state and final_state must have the same number of anyons per qudit"

        nb_qudits = initial_state.nb_qudits
        nb_anyons_per_qudit = initial_state.nb_anyons_per_qudit

        # Get initial and final inputs
        initial_inputs = initial_state.get_inputs()
        final_inputs = final_state.get_inputs()

        # Permute the anyons i and i+1 in the initial inputs
        temp = initial_inputs[braid_index]
        initial_inputs[braid_index] = initial_inputs[braid_index - 1]
        initial_inputs[braid_index - 1] = temp
        # Check if the permuted initial inputs match the final inputs
        if not np.array_equal(initial_inputs, final_inputs):
            print("inputs do not match")
            return 0.0 + 0.0j

        initial_outcomes = deepcopy(
            initial_state.charges[
                nb_qudits * nb_anyons_per_qudit
                + nb_qudits * (nb_anyons_per_qudit - 1) : :
            ]
        )
        final_outcomes = deepcopy(
            final_state.charges[
                nb_qudits * nb_anyons_per_qudit
                + nb_qudits * (nb_anyons_per_qudit - 1) : :
            ]
        )
        if not np.array_equal(initial_outcomes, final_outcomes):
            print("outcomes do not match")
            return 0.0 + 0.0j

        remainder = braid_index % nb_anyons_per_qudit
        if remainder > 0:
            # Braiding within a qudit
            print("braiding within a qudit")
            qudit_index = braid_index // nb_anyons_per_qudit
            """
            [ B^{a(m,r-1), a(m,r)}_{i(m,r-1), a(m,r), a(m,r+1)} ]^{i(m,r)}_{i'(m,r)}
            r = remainder
            """
            for qudit in range(nb_qudits):
                if not np.array_equal(
                    initial_state.get_qudit_state(qudit),
                    final_state.get_qudit_state(qudit),
                ):
                    if qudit != qudit_index:
                        print(initial_state.get_qudit_state(qudit).inputs, 
                              final_state.get_qudit_state(qudit).inputs)
                        print(initial_state.get_qudit_state(qudit).outcomes,
                              final_state.get_qudit_state(qudit).outcomes)
                        print(f"qudit state {qudit} does not much i != f")
                        return 0.0 + 0.0j

            # create standard basis states for initial and final single qudit states
            amplitude = self.compute_standard_braid_component(
                initial_state.get_qudit_state(qudit_index),
                remainder,
                final_state.get_qudit_state(qudit_index),
            )
            return amplitude
        else:
            # Braiding between two qudits
            print("braiding between two qudits")
            first_qudit_index = (braid_index // nb_anyons_per_qudit) - 1
            second_qudit_index = braid_index // nb_anyons_per_qudit
            print(f"first_qudit_index: {first_qudit_index}")

            """
            [K_{
            a(m,q) a(m+1, 0), ..., a(m+1, q), 
            i(m,q-1)}^{j(m-2), j(m)
            }]^{
            j(m-1) i(m,q)i(m+1,1) ... i(m+1,q)
            }_{
            j'(m-1), i'(m,q), i'(m+1,1) ... i'(m+1,q)
            }

            q = nb_anyons_per_qudit - 1
            """
            q_ = nb_anyons_per_qudit - 1
            # a charges of the state a(m,q) a(m+1, 0), ..., a(m+1, q),
            a = []
            m = first_qudit_index
            a.append(deepcopy(initial_state.charges[(m + 1) * q_ - 1]))
            for r in range(0, q_ + 1):
                a.append(deepcopy(initial_state.charges[(m + 1) * q_ + r]))
            # i charges of the state i(m,q-1), i(m,q)i(m+1,0) ... i(m+1,q)
            i = []
            # i(m,q-1)
            i.append(deepcopy(
                initial_state.charges[
                    nb_qudits * nb_anyons_per_qudit
                    + (m + 1) * (nb_anyons_per_qudit - 1)
                    - 2
                ]
            ))
            # i(m,q)
            i.append(
                deepcopy(initial_state.charges[
                    nb_qudits * nb_anyons_per_qudit
                    + (m + 1) * (nb_anyons_per_qudit - 1)
                    - 1
                ])
            )
            # i(m+1,0) ... i(m+1,q)
            for r in range(1, q_ + 1):
                i.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + (m + 1) * (nb_anyons_per_qudit - 1)
                        + r
                    ])
                )

            # i_prime charges of the final state i'(m,q),i'(m+1,0) ... i'(m+1,q)
            i_prime = []
            # i'(m,q)
            i_prime.append(
                deepcopy(final_state.charges[
                    nb_qudits * nb_anyons_per_qudit
                    + (m + 1) * (nb_anyons_per_qudit - 1)
                    - 1
                ])
            )
            # i'(m+1,0) ... i'(m+1,q)
            for r in range(1, q_ + 1):
                i_prime.append(
                    deepcopy(final_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + (m + 1) * (nb_anyons_per_qudit - 1)
                        + r
                    ])
                )
            # root j charges of the state
            # j(m-2), j(m-1), j(m)
            # j are indiced from 0 to nb_qudits - 1
            # (nb_qudits - 1) j indices total
            j = []
            if m == 0:
                j.append(0)  # dummy value for j(m-2)
                j.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_anyons_per_qudit
                        - 1
                        - 1
                    ])
                )  # j(m - 1)
                j.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_qudits * (nb_anyons_per_qudit - 1)
                        + 0
                    ])
                )  # j(m)
            elif m == 1:
                j.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_anyons_per_qudit
                        - 1
                        - 1
                    ])
                )  # j(m - 2)
                j.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_qudits * (nb_anyons_per_qudit - 1)
                        + 0
                    ])
                )  # j(m - 1)
                j.append(
                    deepcopy(initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_qudits * (nb_anyons_per_qudit - 1)
                        + 1
                    ])
                )  # j(m)
            else:
                for r in [m - 2, m - 1, m]:
                    j.append(
                        deepcopy(initial_state.charges[
                            nb_qudits * nb_anyons_per_qudit
                            + nb_qudits * (nb_anyons_per_qudit - 1)
                            + r
                        ])
                    )

            # j'(m-1)
            # j' are indiced from 1 to nb_qudits - 1
            j_prime = []
            if m == 0:
                j_prime.append(
                    deepcopy(final_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_anyons_per_qudit
                        - 1
                        - 1
                    ])
                )
            else:
                j_prime.append(
                    deepcopy(final_state.charges[
                        nb_qudits * nb_anyons_per_qudit
                        + nb_qudits * (nb_anyons_per_qudit - 1)
                        + (m - 1)
                        - 1
                    ])
                )

            """
            Return
            [K_{
            a(m,q) a(m+1, 0), ..., a(m+1, q), 
            i(m,q-1)}^{j(m-2), j(m)
            }]^{
            j(m-1) i(m,q)i(m+1,1) ... i(m+1,q)
            }_{
            j'(m-1), i'(m,q), i'(m+1,1) ... i'(m+1,q)
            }
            """
            knitting_matrix = self.compute_knitting_matrix(
                q=nb_anyons_per_qudit - 1
            )
            # print(knitting_matrix)

            """
            f"a(m,{q})",
            *[f"a(m+1,{r})" for r in range(0, q + 1)],
            f"i(m,{q-1})",
            "j(m-2)",
            "j(m)",
            "j(m-1)",
            f"i(m,{q})",
            *[f"i(m+1,{r})" for r in range(1, q + 1)],
            "jp(m-1)",
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(1, q + 1)],
            """

            return knitting_matrix[
                *a, i[0], j[0], j[2], j[1], *i[1::], j_prime[0], *i_prime
            ]
    
    def generate_computational_braiding_operator(self, index: int, basis: List[AnyonState]):
        """Generates the braiding operator of index 'index' for a system of
        a given number of qudits and anyons per qudit.
        This operator braids anyons at positions 'index' and 'index'+1.

        Parameters
        ----------
        index : int
            The operator's index.
        basis : list of AnyonState objects
            The basis in question.

        Returns
        -------
        List
            Matrix representation of the braiding operator.

        """
        dim = len(basis)
        nb_anyons = basis[0].nb_qudits * basis[0].nb_anyons_per_qudit
        sigmas = np.zeros((dim, dim), dtype=np.complex128)
        for f, base_f in enumerate(basis):
            for i, base_i in enumerate(basis):
                sigmas[f, i] = self.compute_sparse_braid_inner_product(
                    SparseAnyonState(
                        np.concatenate((np.ones(nb_anyons, dtype=int) * base_i.input_charge, base_i.charges)), 
                        base_i.nb_qudits,
                        base_i.nb_anyons_per_qudit
                    ), index, 
                    SparseAnyonState(
                        np.concatenate((np.ones(nb_anyons, dtype=int) * base_f.input_charge, base_f.charges)), 
                        base_f.nb_qudits,
                        base_f.nb_anyons_per_qudit
                    )
                    )

        return sigmas
