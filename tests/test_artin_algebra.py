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
import pytest

from tqsim import AnyonicCircuit
from tqsim.models.fibonacci import FIBONACCI_MODEL
from tqsim.models.ising import ISING_MODEL


class TestFibonacciArtinAlgebra:
    """Test Artin algebra relations for Fibonacci anyons."""

    @pytest.mark.parametrize("nb_qudits", [1, 2])
    def test_commutation_relations(self, nb_qudits: int) -> None:
        r"""Test that braid operators commute when |i-j| > 1.

        For braid indices i and j where |i-j| > 1, the braid operators
        should satisfy: σ_i σ_j = σ_j σ_i
        """
        nb_anyons_per_qudit = 3
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=FIBONACCI_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test commutation for all pairs where |i-j| > 1
        for braid_1 in range(nb_braid_generators):
            for braid_2 in range(braid_1 + 2, nb_braid_generators):
                left = (
                    circuit.braiding_operators[braid_1]
                    @ circuit.braiding_operators[braid_2]
                )
                right = (
                    circuit.braiding_operators[braid_2]
                    @ circuit.braiding_operators[braid_1]
                )
                diff_norm = np.linalg.norm(left - right)
                assert np.isclose(
                    diff_norm, 0, atol=1e-14
                ), f"Commutation failed for braids {braid_1} and {braid_2}"

    @pytest.mark.parametrize("nb_qudits", [1, 2])
    def test_yang_baxter_relation(self, nb_qudits: int) -> None:
        r"""Test the Yang-Baxter (Artin braid) relation.

        For adjacent braid indices i and i+1, the operators should satisfy:
        σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}
        """
        nb_anyons_per_qudit = 3
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=FIBONACCI_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test Yang-Baxter relation for all adjacent pairs
        for braid in range(nb_braid_generators - 1):
            left = (
                circuit.braiding_operators[braid]
                @ circuit.braiding_operators[braid + 1]
                @ circuit.braiding_operators[braid]
            )
            right = (
                circuit.braiding_operators[braid + 1]
                @ circuit.braiding_operators[braid]
                @ circuit.braiding_operators[braid + 1]
            )
            diff_norm = np.linalg.norm(left - right)
            assert np.isclose(
                diff_norm, 0, atol=1e-14
            ), f"Yang-Baxter relation failed for braid {braid}"

    @pytest.mark.parametrize("nb_qudits", [1, 2])
    def test_unitarity(self, nb_qudits: int) -> None:
        r"""Test that braid operators are unitary.

        Each braid operator σ_i should satisfy:
        σ_i σ_i^† = σ_i^† σ_i = I
        """
        nb_anyons_per_qudit = 3
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=FIBONACCI_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test unitarity for all braid operators
        for braid in range(nb_braid_generators):
            operator = circuit.braiding_operators[braid]
            inverse = np.linalg.inv(operator)
            identity = np.eye(operator.shape[0])

            # Test σ_i σ_i^{-1} = I
            left = operator @ inverse
            assert np.isclose(
                np.linalg.norm(left - identity), 0, atol=1e-14
            ), f"Left unitarity failed for braid {braid}"

            # Test σ_i^{-1} σ_i = I
            right = inverse @ operator
            assert np.isclose(
                np.linalg.norm(right - identity), 0, atol=1e-14
            ), f"Right unitarity failed for braid {braid}"

            # Test left = right (commutativity of operator with its inverse)
            assert np.isclose(
                np.linalg.norm(left - right), 0, atol=1e-14
            ), f"Inverse commutativity failed for braid {braid}"


class TestIsingArtinAlgebra:
    """Test Artin algebra relations for Ising anyons."""

    @pytest.mark.parametrize("nb_anyons_per_qudit", [3, 4, 5])
    def test_commutation_relations(self, nb_anyons_per_qudit: int) -> None:
        r"""Test that braid operators commute when |i-j| > 1.

        For braid indices i and j where |i-j| > 1, the braid operators
        should satisfy: σ_i σ_j = σ_j σ_i
        """
        nb_qudits = 1
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=ISING_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test commutation for all pairs where |i-j| > 1
        for braid_1 in range(nb_braid_generators):
            for braid_2 in range(braid_1 + 2, nb_braid_generators):
                left = (
                    circuit.braiding_operators[braid_1]
                    @ circuit.braiding_operators[braid_2]
                )
                right = (
                    circuit.braiding_operators[braid_2]
                    @ circuit.braiding_operators[braid_1]
                )
                diff_norm = np.linalg.norm(left - right)
                assert np.isclose(
                    diff_norm, 0, atol=1e-14
                ), f"Commutation failed for braids {braid_1} and {braid_2}"

    @pytest.mark.parametrize("nb_anyons_per_qudit", [3, 4, 5])
    def test_yang_baxter_relation(self, nb_anyons_per_qudit: int) -> None:
        r"""Test the Yang-Baxter (Artin braid) relation.

        For adjacent braid indices i and i+1, the operators should satisfy:
        σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}
        """
        nb_qudits = 1
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=ISING_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test Yang-Baxter relation for all adjacent pairs
        for braid in range(nb_braid_generators - 1):
            left = (
                circuit.braiding_operators[braid]
                @ circuit.braiding_operators[braid + 1]
                @ circuit.braiding_operators[braid]
            )
            right = (
                circuit.braiding_operators[braid + 1]
                @ circuit.braiding_operators[braid]
                @ circuit.braiding_operators[braid + 1]
            )
            diff_norm = np.linalg.norm(left - right)
            assert np.isclose(
                diff_norm, 0, atol=1e-14
            ), f"Yang-Baxter relation failed for braid {braid}"

    @pytest.mark.parametrize("nb_anyons_per_qudit", [3, 4, 5])
    def test_unitarity(self, nb_anyons_per_qudit: int) -> None:
        r"""Test that braid operators are unitary.

        Each braid operator σ_i should satisfy:
        σ_i σ_i^† = σ_i^† σ_i = I
        """
        nb_qudits = 1
        circuit = AnyonicCircuit(
            nb_qudits, nb_anyons_per_qudit, model=ISING_MODEL, input_charge=1
        )

        nb_anyons = nb_qudits * nb_anyons_per_qudit
        nb_braid_generators = nb_anyons - 1

        # Test unitarity for all braid operators
        for braid in range(nb_braid_generators):
            operator = circuit.braiding_operators[braid]
            inverse = np.linalg.inv(operator)
            identity = np.eye(operator.shape[0])

            # Test σ_i σ_i^{-1} = I
            left = operator @ inverse
            assert np.isclose(
                np.linalg.norm(left - identity), 0, atol=1e-14
            ), f"Left unitarity failed for braid {braid}"

            # Test σ_i^{-1} σ_i = I
            right = inverse @ operator
            assert np.isclose(
                np.linalg.norm(right - identity), 0, atol=1e-14
            ), f"Right unitarity failed for braid {braid}"

            # Test left = right (commutativity of operator with its inverse)
            assert np.isclose(
                np.linalg.norm(left - right), 0, atol=1e-14
            ), f"Inverse commutativity failed for braid {braid}"
