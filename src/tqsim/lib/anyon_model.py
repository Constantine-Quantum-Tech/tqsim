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

import os
from copy import deepcopy

import numpy as np

from tqsim.config import STORE_PATH
from tqsim.lib.anyon_state import (
    AnyonState,
    SparseAnyonState,
    StandardAnyonState,
)
from tqsim.lib.utils import einsum_with_names


class AnyonModel:
    """Class representing an anyon model.

    Parameters
    ----------
    N_symbols : np.ndarray
        The fusion rules N symbols tensor.
    F_matrix : np.ndarray
        The F matrix tensor.
    R_matrix : np.ndarray
        The R matrix tensor.
    name : str, optional
        The name of the anyon model. If None, a random name is generated.
    force_recache : bool, optional
        If True, forces the regeneration of cached matrices. Default is False.

    Attributes
    ----------
    N_symbols : np.ndarray
        The fusion rules N symbols tensor.
    F_matrix : np.ndarray
        The F matrix tensor.
    R_matrix : np.ndarray
        The R matrix tensor.
    name : str
        The name of the anyon model.

    Example 1:
        >>> from tqsim.models.fibonacci import FIBONACCI_MODEL
        >>> model = FIBONACCI_MODEL
        >>> print(model.n_symbols)
        >>> print(model.f_matrix)
        >>> print(model.r_matrix)

    Example 2: Z_N model (Abelian model)
        >>> N = 5
        >>> n_symbols = np.zeros((N, N, N), dtype=int)
        >>> for i, j, k in itertools.product(range(N), repeat=3):
        >>>     if (i + j) % N == k:
        >>>         n_symbols[i, j, k] = 1

        >>> f_matrix = np.zeros((N, N, N, N, N, N), dtype=complex)
        >>> for i, j, k, fusion_result, m, n in itertools.product(range(N), repeat=6):
        >>>     if (i + j + k) % N == fusion_result and (i + j) % N == m and (j + k) % N == n:
        >>>         f_matrix[i, j, k, fusion_result, m, n] = 1

        >>> r_matrix = np.zeros((N, N, N), dtype=complex)
        >>> for i, j, k in itertools.product(range(N), repeat=3):
        >>>     if (i + j) % N == k:
        >>>         r_matrix[i, j, k] = np.exp(2j * np.pi * i * j / N)

        >>> zn_model = AnyonModel(n_symbols, f_matrix, r_matrix, name="Z_N")
    """

    def __init__(
        self,
        n_symbols: np.ndarray,
        f_matrix: np.ndarray,
        r_matrix: np.ndarray,
        name=None,
        force_recache=False,
    ):
        assert n_symbols.ndim == 3, "N_symbols must be a 3D tensor"
        assert f_matrix.ndim == 6, "F_matrix must be a 6D tensor"
        assert r_matrix.ndim == 3, "R_matrix must be a 5D tensor"

        if force_recache:
            folder_path = os.path.join(STORE_PATH, f"{name}-q-*")
            if os.path.exists(folder_path):
                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)

        self._n_symbols = n_symbols
        self._f_matrix = f_matrix
        self._r_matrix = r_matrix
        self._b_matrix = self._compute_braiding_matrix()
        self._k_matrices = {}
        self.nb_charges = self.n_symbols.shape[0]

        if name is None:
            self._name = f"model-{np.random.randint(1000)}"
        else:
            self._name = name

    @property
    def n_symbols(self):
        """Returns the fusion rules N symbols tensor.

        Returns
        -------
        np.ndarray
            The N symbols tensor.

        """
        return self._n_symbols

    @property
    def f_matrix(self):
        """Returns the F matrix tensor.

        Returns
        -------
        np.ndarray
            The F matrix tensor.

        """
        return self._f_matrix

    @property
    def r_matrix(self):
        """Returns the R matrix tensor.

        Returns
        -------
        np.ndarray
            The R matrix tensor.

        """
        return self._r_matrix

    @property
    def b_matrix(self):
        """Returns the braiding matrix tensor.

        Returns
        -------
        np.ndarray
            The braiding matrix tensor.

        return self._b_matrix
            The braiding matrix tensor.

        """
        return self._b_matrix

    @property
    def k_matrices(self):
        """Returns the K matrices dictionary.

        Returns
        -------
        dict
            The K matrices dictionary.

        """
        return self._k_matrices

    @property
    def name(self):
        """Returns the name of the anyon model.

        Returns
        -------
        str
            The name of the anyon model.

        """
        return self._name

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
        # check that n_symbols[anyon1[i], anyon2[i], outcome[i]] == 1 for all i
        return self.n_symbols[anyon1, anyon2, outcome] == np.ones_like(anyon1)

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
            self.f_matrix,
            self.r_matrix,
            self.f_matrix.conjugate(),
        )
        return b_matrix

    def _compute_l_matrix(self, q: int):
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
            "q must be strictly positive. " "For q=1, L is just the braiding matrix."
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
            f_dag = np.conjugate(self.f_matrix).swapaxes(-1, -2)
            terms.append((f_dag, labels))
            r += 1

        # --- B tensor
        b_labels = (
            f"i(m,{q-1})",
            f"a(m,{q})",
            "a(m+1,0)",
            "p(1)",
            f"i(m,{q})",
            f"ip(m,{q})",
        )
        terms.append((self.b_matrix, b_labels))

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
            terms.append((self.f_matrix, labels))
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

    def compute_knitting_matrix(self, q: int, return_l=False):
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
        # Check if K matrix is already computed
        folder_path = os.path.join(STORE_PATH, f"{self.name}-q-{q}")
        filename = os.path.join(folder_path, "-knitting-matrix.npy")
        if os.path.exists(filename):
            return np.load(filename)

        # Compute K matrix
        assert q > 0, (
            "q must be strictly positive. " "For q=0, K is just the braiding matrix."
        )

        terms = []

        # --- Left F factor ---
        # Example label order (must match how F is actually stored):
        f_labels = (
            "j(m-2)",
            f"i(m,{q})",
            f"i(m+1,{q})",
            "j(m)",
            "j(m-1)",
            "k",
        )
        terms.append((self.f_matrix, f_labels))

        # --- L factor ---
        # L already carries many indices (a, i, i', ...).
        # Here we assume l_labels is known / fixed.
        # Example layout (you must adapt to your actual storage order!):
        # out_labels = (
        #     [f"a(m,{q})"] + [f"a(m+1,{r})" for r in range(0, q+1)] +
        #     [f"i(m,{q-1})"] + [f"p({q+1})"] +
        #     [f"i(m,{q})"] + [f"i(m+1,{r})" for r in range(0, q+1)] +
        #     [f"ip(m,{q})"] + [f"ip(m+1,{r})" for r in range(0, q+1)]
        # )
        l_labels = (
            f"a(m,{q})",
            *[f"a(m+1,{r})" for r in range(0, q + 1)],
            f"i(m,{q-1})",
            "k",
            f"i(m,{q})",
            *[f"i(m+1,{r})" for r in range(1, q + 1)],
            f"ip(m,{q})",
            *[f"ip(m+1,{r})" for r in range(1, q + 1)],
        )
        l_matrix = self._compute_l_matrix(q)
        terms.append((l_matrix, l_labels))

        # --- Right F dagger factor ---
        # dagger = conjugate and swapaxes(5,4): swap 'k' and 'j(m-1)'
        f_dag = np.conjugate(np.swapaxes(self.f_matrix, 5, 4))
        f_labels_dag = (
            "j(m-2)",
            f"ip(m,{q})",
            f"ip(m+1,{q})",
            "j(m)",
            "k",
            "jp(m-1)",
        )
        terms.append((f_dag, f_labels_dag))

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
        k = einsum_with_names(terms, out_labels)
        self._k_matrices[q] = k

        # Store K matrix to file
        os.makedirs(folder_path, exist_ok=True)
        np.save(filename, k)

        if return_l:
            return k, l_matrix
        return k

    def check_state(self, state) -> bool:
        if isinstance(state, SparseAnyonState):
            return self.check_sparse_state(
                state.charges, state.nb_qudits, state.nb_anyons_per_qudit
            )
        elif isinstance(state, StandardAnyonState):
            return self.check_standard_basis_state(state.inputs, state.outcomes)
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
            start = nb_qudits * nb_anyons_per_qudit + i * (nb_anyons_per_qudit - 1)
            end = nb_qudits * nb_anyons_per_qudit + (i + 1) * (nb_anyons_per_qudit - 1)
            outcomes = charges[start:end]
            check = self.check_standard_basis_state(inputs, outcomes)
            if not check:
                return False

        if nb_qudits == 1:
            return True

        indices = [
            nb_qudits * nb_anyons_per_qudit + (i + 1) * (nb_anyons_per_qudit - 1) - 1
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

        initial_inputs = deepcopy(initial_state.inputs)
        final_inputs = deepcopy(final_state.inputs)

        # Permute the anyons i and i+1 in the initial inputs
        temp = deepcopy(initial_inputs[braid_index])
        initial_inputs[braid_index] = deepcopy(initial_inputs[braid_index - 1])
        initial_inputs[braid_index - 1] = deepcopy(temp)
        # Check if the permuted initial inputs match the final inputs
        if not np.array_equal(initial_inputs, final_inputs):
            return 0.0 + 0.0j

        initial_outcomes = deepcopy(initial_state.outcomes)
        final_outcomes = deepcopy(final_state.outcomes)
        if braid_index >= 2:
            final_outcomes[braid_index - 2] = deepcopy(
                initial_outcomes[braid_index - 2]
            )

        if not np.array_equal(initial_outcomes, final_outcomes):
            return 0.0 + 0.0j

        if braid_index == 1:
            # Braiding the first two anyons
            a = 0  # vacuum
        elif braid_index == 2:
            a = deepcopy(initial_state.inputs[0])
        else:
            a = deepcopy(initial_state.outcomes[braid_index - 3])

        b = deepcopy(initial_state.inputs[braid_index - 1])

        c = deepcopy(initial_state.inputs[braid_index])

        if braid_index == 1:
            i = deepcopy(initial_state.inputs[0])
            m = deepcopy(final_state.inputs[0])
        else:
            i = deepcopy(initial_state.outcomes[braid_index - 2])
            m = deepcopy(final_state.outcomes[braid_index - 2])

        j = deepcopy(final_state.outcomes[braid_index - 1])

        amplitude = self.b_matrix[a, b, c, j, i, m]
        return amplitude

    def _validate_sparse_braid_inputs(
        self,
        initial_state: SparseAnyonState,
        braid_index: int,
        final_state: SparseAnyonState,
    ):
        """Validate inputs for sparse braid inner product computation."""
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
            initial_state.nb_anyons_per_qudit == final_state.nb_anyons_per_qudit
        ), "initial_state and final_state must have the same number of anyons per qudit"

    def _check_permuted_inputs_match(
        self,
        initial_state: SparseAnyonState,
        final_state: SparseAnyonState,
        braid_index: int,
    ):
        """Check if permuted initial inputs match final inputs."""
        initial_inputs = deepcopy(initial_state.get_inputs())
        final_inputs = deepcopy(final_state.get_inputs())

        # Permute the anyons i and i+1 in the initial inputs
        temp = deepcopy(initial_inputs[braid_index])
        initial_inputs[braid_index] = deepcopy(initial_inputs[braid_index - 1])
        initial_inputs[braid_index - 1] = deepcopy(temp)

        return np.array_equal(initial_inputs, final_inputs)

    def _compute_within_qudit_braid(
        self,
        initial_state: SparseAnyonState,
        final_state: SparseAnyonState,
        qudit_index: int,
        remainder: int,
        nb_qudits: int,
    ):
        """Compute braiding within a single qudit."""
        for qudit in range(nb_qudits):
            if not np.array_equal(
                initial_state.get_qudit_state(qudit),
                final_state.get_qudit_state(qudit),
            ):
                if qudit != qudit_index:
                    return 0.0 + 0.0j

        initial_outcomes = deepcopy(initial_state.get_outcomes())
        final_outcomes = deepcopy(final_state.get_outcomes())
        if not np.array_equal(initial_outcomes, final_outcomes):
            return 0.0 + 0.0j

        # Check that constant nodes stay fixed
        qubit_state_initial = deepcopy(initial_state.get_qudit_state(qudit_index))
        qubit_state_final = deepcopy(final_state.get_qudit_state(qudit_index))

        # create standard basis states for initial and final single qudit states
        return self.compute_standard_braid_component(
            qubit_state_initial,
            remainder,
            qubit_state_final,
        )

    def _validate_qudit_states_between_qudits(
        self,
        initial_state: SparseAnyonState,
        final_state: SparseAnyonState,
        first_qudit_index: int,
        second_qudit_index: int,
        nb_qudits: int,
    ):
        """Validate that only the two braided qudits can differ."""
        for qudit in range(nb_qudits):
            if not np.array_equal(
                initial_state.get_qudit_state(qudit),
                final_state.get_qudit_state(qudit),
            ):
                if qudit not in [first_qudit_index, second_qudit_index]:
                    return False
        return True

    def _validate_outcomes_between_qudits(
        self,
        initial_state: SparseAnyonState,
        final_state: SparseAnyonState,
        m: int,
        nb_qudits: int,
        nb_anyons_per_qudit: int,
    ):
        """Validate outcomes for braiding between qudits."""
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

        if m > 0:
            final_outcomes[m - 1] = initial_outcomes[m - 1]

        if not np.array_equal(initial_outcomes, final_outcomes):
            return False

        unmodified_i_initial = initial_state.charges[
            nb_qudits * nb_anyons_per_qudit
            + m * (nb_anyons_per_qudit - 1) : nb_qudits * nb_anyons_per_qudit
            + m * (nb_anyons_per_qudit - 1)
            + nb_anyons_per_qudit
            - 2
        ]

        unmodified_i_final = final_state.charges[
            nb_qudits * nb_anyons_per_qudit
            + m * (nb_anyons_per_qudit - 1) : nb_qudits * nb_anyons_per_qudit
            + m * (nb_anyons_per_qudit - 1)
            + nb_anyons_per_qudit
            - 2
        ]

        return np.array_equal(unmodified_i_initial, unmodified_i_final)

    def _extract_j_charges(
        self,
        initial_state: SparseAnyonState,
        m: int,
        q: int,
        nb_qudits: int,
        nb_anyons_per_qudit: int,
    ):
        """Extract j charges for knitting matrix computation."""
        j = []
        if m == 0:
            j.append(0)  # dummy value for j(m-2)
            j.append(
                deepcopy(initial_state.charges[nb_qudits * nb_anyons_per_qudit + q - 1])
            )  # j(m - 1)
            j.append(
                deepcopy(
                    initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit + nb_qudits * q
                    ]
                )
            )  # j(m)
        elif m == 1:
            j.append(
                deepcopy(initial_state.charges[nb_qudits * nb_anyons_per_qudit + q - 1])
            )  # j(m - 2)
            j.append(
                deepcopy(
                    initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit + nb_qudits * q
                    ]
                )
            )  # j(m - 1)
            j.append(
                deepcopy(
                    initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit + nb_qudits * q + 1
                    ]
                )
            )  # j(m)
        else:
            for r in [m - 2, m - 1, m]:
                j.append(
                    deepcopy(
                        initial_state.charges[
                            nb_qudits * nb_anyons_per_qudit + nb_qudits * q + r
                        ]
                    )
                )
        return j

    def _compute_between_qudits_braid(
        self,
        initial_state: SparseAnyonState,
        final_state: SparseAnyonState,
        first_qudit_index: int,
        second_qudit_index: int,
        nb_qudits: int,
        nb_anyons_per_qudit: int,
    ):
        """Compute braiding between two qudits."""
        m = first_qudit_index

        if not self._validate_qudit_states_between_qudits(
            initial_state, final_state, first_qudit_index, second_qudit_index, nb_qudits
        ):
            return 0.0 + 0.0j

        if not self._validate_outcomes_between_qudits(
            initial_state, final_state, m, nb_qudits, nb_anyons_per_qudit
        ):
            return 0.0 + 0.0j

        q = nb_anyons_per_qudit - 1

        # a charges of the state a(m,q) a(m+1, 0), ..., a(m+1, q),
        a = []
        # a(m,q)
        a.append(deepcopy(initial_state.charges[(m + 1) * nb_anyons_per_qudit - 1]))
        # a(m+1,0) ... a(m+1,q)
        for r in range(0, q + 1):
            a.append(deepcopy(initial_state.charges[(m + 1) * nb_anyons_per_qudit + r]))

        # i charges of the state i(m,q-1), i(m,q)i(m+1,1) ... i(m+1,q)
        i = []
        # i(m,q-1)
        i.append(
            deepcopy(
                initial_state.charges[nb_qudits * nb_anyons_per_qudit + (m + 1) * q - 2]
            )
        )
        # i(m,q)
        i.append(
            deepcopy(
                initial_state.charges[nb_qudits * nb_anyons_per_qudit + (m + 1) * q - 1]
            )
        )
        # i(m+1,1) ... i(m+1,q)
        for r in range(1, q + 1):
            i.append(
                deepcopy(
                    initial_state.charges[
                        nb_qudits * nb_anyons_per_qudit + (m + 1) * q + r - 1
                    ]
                )
            )

        # i_prime charges of the final state i'(m,q),i'(m+1,0) ... i'(m+1,q)
        i_prime = []
        # i'(m,q)
        i_prime.append(
            deepcopy(
                final_state.charges[nb_qudits * nb_anyons_per_qudit + (m + 1) * q - 1]
            )
        )
        # i'(m+1,1) ... i'(m+1,q)
        for r in range(1, q + 1):
            i_prime.append(
                deepcopy(
                    final_state.charges[
                        nb_qudits * nb_anyons_per_qudit + (m + 1) * q + r - 1
                    ]
                )
            )

        # root j charges of the state
        j = self._extract_j_charges(initial_state, m, q, nb_qudits, nb_anyons_per_qudit)

        # j'(m-1)
        j_prime = []
        if m == 0:
            j_prime.append(
                deepcopy(final_state.charges[nb_qudits * nb_anyons_per_qudit + q - 1])
            )
        else:
            j_prime.append(
                deepcopy(
                    final_state.charges[
                        nb_qudits * nb_anyons_per_qudit + nb_qudits * q + (m - 1)
                    ]
                )
            )

        knitting_matrix = self._k_matrices.get(q, self.compute_knitting_matrix(q=q))

        return knitting_matrix[
            *a, i[0], j[0], j[2], j[1], *i[1::], j_prime[0], *i_prime
        ]

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
        self._validate_sparse_braid_inputs(initial_state, braid_index, final_state)

        nb_qudits = initial_state.nb_qudits
        nb_anyons_per_qudit = initial_state.nb_anyons_per_qudit

        # Check if the permuted initial inputs match the final inputs
        if not self._check_permuted_inputs_match(
            initial_state, final_state, braid_index
        ):
            return 0.0 + 0.0j

        remainder = braid_index % nb_anyons_per_qudit
        if remainder > 0:
            # Braiding within a qudit
            qudit_index = braid_index // nb_anyons_per_qudit
            return self._compute_within_qudit_braid(
                initial_state, final_state, qudit_index, remainder, nb_qudits
            )
        else:
            # Braiding between two qudits
            first_qudit_index = (braid_index // nb_anyons_per_qudit) - 1
            second_qudit_index = braid_index // nb_anyons_per_qudit
            return self._compute_between_qudits_braid(
                initial_state,
                final_state,
                first_qudit_index,
                second_qudit_index,
                nb_qudits,
                nb_anyons_per_qudit,
            )

    def generate_computational_braiding_operator(
        self, index: int, basis: list[AnyonState]
    ):
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
                        np.concatenate(
                            (
                                np.ones(nb_anyons, dtype=int) * base_i.input_charge,
                                base_i.charges,
                            )
                        ),
                        base_i.nb_qudits,
                        base_i.nb_anyons_per_qudit,
                    ),
                    index,
                    SparseAnyonState(
                        np.concatenate(
                            (
                                np.ones(nb_anyons, dtype=int) * base_f.input_charge,
                                base_f.charges,
                            )
                        ),
                        base_f.nb_qudits,
                        base_f.nb_anyons_per_qudit,
                    ),
                )

        return sigmas
