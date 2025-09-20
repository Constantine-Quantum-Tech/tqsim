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

import numpy as np
from typing import Tuple
from lib.utils import einsum_with_names
from lib.anyon_state import StandardAnyonState, SparseAnyonState


class AnyonModel:

    def __init__(self, N_symbols, F_matrix, R_matrix):
        self.N_symbols = N_symbols
        self.F_matrix = F_matrix
        self.R_matrix = R_matrix
        self.braiding_matrix = self._compute_braiding_matrix()

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
        [L_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{p(q+1)}]^{
        i(m,q) i(m+1,0) ... i(m+1,q)
        }_{
        i'(m,q), i'(m+1,0) ... i'(m+1,q)
        }
        =
        \sum_{p(1), .., p(q), i'(m,q), i'(m+1,1), .., i'(m+1,q)}
        prod_{r=1}^{q} [ F_{i(m, q), i(m+1, q-r),
                         a(m+1, q-r+1)}^{p(q-r+2)}
                        ].dag^{i(m+1, q-r+1)}_{p(q-r+1)}
        [ B^{p(1)}_{i(m,q-1), a(m, q), a(m+1, 0)} ]^{i(m, q)}_{i'(m, q)}
        prod_{r=1}^{q} [ F_{i'(m, q), i'(m+1, q-r),
                         a(m+1, q-r+1)}^{p(q-r+2)}
                        ]^{p(q-r+1)}_{i'(m+1, q-r+1)}
        """
        assert (
            q > 0
        ), ("q must be strictly positive. "
            "For q=0, L is just the braiding matrix.")
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
            # conj + swap last two axes for dagger
            F_dag = np.conjugate(self.F_matrix).swapaxes(-1, -2)
            terms.append((F_dag, labels))
            r += 1

        # --- B tensor
        B_labels = (
            "p(1)",
            f"i(m,{q-1})",
            f"a(m,{q})",
            "a(m+1,0)",
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
            + [f"i(m+1,{r})" for r in range(0, q + 1)]
            + [f"ip(m,{q})"]
            + [f"ip(m+1,{r})" for r in range(0, q + 1)]
        )

        # Call the helper from earlier
        return einsum_with_names(terms, out_labels)

    def compute_knitting_matrix(self, q: int):
        r"""
        See Appendix of https://arxiv.org/abs/2307.01892

        [K_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{j(m-2), j(m)}]^{
        j(m-1) i(m,q)i(m+1,0) ... i(m+1,q)
        }_{
        j'(m-1), i'(m,q), i'(m+1,0) ... i'(m+1,q)
        }
        =
        \sum_{k}
        [ F^{j(m)}_{j(m-2), i(m,q), i(m+1,q)} ]^{j(m-1)}_{k}
        [L_{a(m,q) a(m+1, 0), ..., a(m+1, q), i(m,q-1)}^{k}]^{
        i(m,q) i(m+1,0) ... i(m+1,q)
        }_{
        i'(m,q), i'(m+1,0) ... i'(m+1,q)
        }
        [ F^{j(m)}_{j(m-2), i'(m,q), i'(m+1,q)} ].dagger^{k}_{j'(m-1)}

        """
        assert (
            q > 0
        ), ("q must be strictly positive. "
            "For q=0, K is just the braiding matrix.")
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
            *[f"i(m+1,{r})" for r in range(0, q + 1)],
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(0, q + 1)],
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
            *[f"i(m+1,{r})" for r in range(0, q + 1)],
            "jp(m-1)",
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(0, q + 1)],
        )

        # Perform contraction
        K = einsum_with_names(terms, out_labels)
        return K
    
    def check_state(self, state) -> bool:
        if isinstance(state, SparseAnyonState):
            return self.check_sparse_state(state.charges, state.nb_qudits, state.nb_anyons_per_qudit)
        elif isinstance(state, StandardAnyonState):
            return self.check_standard_basis_state(state.inputs, state.outcomes)
        else:
            raise ValueError("State must be either SparseAnyonState or StandardAnyonState")
    
    
    def check_standard_basis_state(self, inputs: np.ndarray, outcomes: np.ndarray) -> bool:
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
        assert inputs.shape[0] == outcomes.shape[0] + 1, "inputs must have one more element than outcomes"
        anyon_1 = inputs[1::]
        anyon_2 = np.zeros(len(anyon_1), dtype=int)
        anyon_2[0] = inputs[0]
        anyon_2[1::] = outcomes[0:-1]
        rules = self.check_rule(anyon_1, anyon_2, outcomes)
        return np.all(rules)

    def check_sparse_state(self, charges: np.ndarray, 
                           nb_qudits: int, nb_anyons_per_qudit: int) -> bool:
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
            start = nb_qudits * nb_anyons_per_qudit + i * (nb_anyons_per_qudit - 1)
            end = nb_qudits * nb_anyons_per_qudit + (i + 1) * (nb_anyons_per_qudit - 1)
            outcomes = charges[start:end]
            check = self.check_standard_basis_state(inputs, outcomes)
            if not check:
                return False
        indices = [nb_qudits * nb_anyons_per_qudit + (i+1) * (nb_anyons_per_qudit - 1) - 1 for i in range(nb_qudits)]
        inputs = charges[indices]
        outcomes = charges[nb_qudits * (2 * nb_anyons_per_qudit -1)::]
        return self.check_standard_basis_state(inputs, outcomes)
    
    def check_general_state(self, inputs: Tuple, outcome: int) -> bool:
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
