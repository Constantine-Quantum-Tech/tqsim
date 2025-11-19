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
import pytest  # type: ignore[import-not-found]

from tqsim.lib.anyon_model import AnyonModel


class TestAnyonModel:
    """Test suite for AnyonModel class based on notebook experiments."""

    @pytest.fixture
    def fibonacci_model(self) -> AnyonModel:
        """Create a Fibonacci anyon model as tested in the notebook."""
        # Fusion matrix setup from notebook
        fusion_matrix = self._build_fusion_matrix()
        f_matrix = self._build_f_matrix()
        r_matrix = self._build_r_matrix()

        return AnyonModel(fusion_matrix, f_matrix, r_matrix)

    def _build_fusion_matrix(self) -> np.ndarray:
        """Build the fusion matrix."""
        fusion_matrix = np.zeros((2, 2, 2))
        for a1 in range(2):
            for a2 in range(2):
                for outcome in range(2):
                    if (a1, a2) == (1, 1):
                        fusion_matrix[a1, a2, outcome] = 1
                    else:
                        if (a1 + a2) == outcome:
                            fusion_matrix[a1, a2, outcome] = 1
        return fusion_matrix

    def _build_f_matrix(self) -> np.ndarray:
        """Build the F matrix."""
        f_matrix = np.zeros((2, 2, 2, 2, 2, 2)) * (1 + 0j)
        for a1 in range(2):
            for a2 in range(2):
                for a3 in range(2):
                    for outcome in range(2):
                        f_matrix[a1, a2, a3, outcome] = self._get_f_matrix(
                            a1, a2, a3, outcome
                        )
        return f_matrix

    def _build_r_matrix(self) -> np.ndarray:
        """Build the R matrix."""
        r_matrix = np.zeros((2, 2, 2)) * (1 + 0j)
        for a1 in range(2):
            for a2 in range(2):
                r_matrix[a1, a2] = self._get_r_matrix(a1, a2).diagonal()
        return r_matrix

    def _get_f_matrix_sum_4(self) -> np.ndarray:
        """Get F matrix when sum = 4."""
        inv_phi = (np.sqrt(5) - 1) / 2
        return np.array([[inv_phi, np.sqrt(inv_phi)], [np.sqrt(inv_phi), -inv_phi]])

    def _get_f_matrix_sum_3(self) -> np.ndarray:
        """Get F matrix when sum = 3."""
        return np.array([[0, 0], [0, 1]])

    def _get_f_matrix_sum_2(self, a1: int, a2: int, a3: int, outcome: int) -> np.ndarray:
        """Get F matrix when sum = 2."""
        if a1 + a2 == 2:
            return np.array([[0, 1], [0, 0]])
        elif a2 + a3 == 2:
            return np.array([[0, 0], [1, 0]])
        elif a1 + a3 == 2:
            return np.array([[0, 0], [0, 1]])
        elif a3 + outcome == 2:
            return np.array([[0, 1], [0, 0]])
        elif a1 + outcome == 2:
            return np.array([[0, 0], [1, 0]])
        elif a2 + outcome == 2:
            return np.array([[0, 0], [0, 1]])
        return np.array([[0, 0], [0, 0]])

    def _get_f_matrix_sum_0(self) -> np.ndarray:
        """Get F matrix when sum = 0."""
        return np.array([[1, 0], [0, 0]])

    def _get_f_matrix(self, a1: int, a2: int, a3: int, outcome: int) -> np.ndarray:
        """F matrix helper function from notebook."""
        total = a1 + a2 + a3 + outcome

        if total == 4:
            return self._get_f_matrix_sum_4()
        elif total == 3:
            return self._get_f_matrix_sum_3()
        elif total == 2:
            return self._get_f_matrix_sum_2(a1, a2, a3, outcome)
        elif total == 0:
            return self._get_f_matrix_sum_0()
        return np.array([[0, 0], [0, 0]])

    def _get_r_matrix(self, a1: int, a2: int) -> np.ndarray:
        """R matrix helper function from notebook."""
        if a1 + a2 == 2:
            matrix = np.array(
                [
                    [np.exp(-4 * np.pi * 1j / 5), 0],
                    [0, np.exp(3 * np.pi * 1j / 5)],
                ]
            )
        else:
            matrix = np.array([[1, 0], [0, 1]])
        return matrix

    def test_model_initialization(self, fibonacci_model: AnyonModel) -> None:
        """Test that the model initializes correctly."""
        assert fibonacci_model is not None
        assert hasattr(fibonacci_model, "n_symbols")
        assert hasattr(fibonacci_model, "f_matrix")
        assert hasattr(fibonacci_model, "r_matrix")
        assert hasattr(fibonacci_model, "b_matrix")

    def test_braiding_matrix_shape(self, fibonacci_model: AnyonModel) -> None:
        """Test that braiding matrix has correct shape as observed in notebook."""
        b_matrix = fibonacci_model.b_matrix
        expected_shape = (2, 2, 2, 2, 2, 2)
        assert (
            b_matrix.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {b_matrix.shape}"

    def test_l_matrix_computation_q1(self, fibonacci_model: AnyonModel) -> None:
        """Test L matrix computation for q=1 (successful case from notebook)."""
        l1 = fibonacci_model._compute_l_matrix(q=1)
        expected_shape = (2, 2, 2, 2, 2, 2, 2, 2, 2)
        assert (
            l1.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {l1.shape}"

    def test_l_matrix_computation_q0_raises_error(self, fibonacci_model: AnyonModel) -> None:
        """Test that L matrix computation for q=0 raises AssertionError as observed in notebook."""
        with pytest.raises(AssertionError, match="q must be strictly positive"):
            fibonacci_model._compute_l_matrix(q=0)

    def test_knitting_matrix_computation_q1(self, fibonacci_model: AnyonModel) -> None:
        """Test knitting matrix computation for q=1 (successful case from notebook)."""
        k1 = fibonacci_model.compute_knitting_matrix(q=1, return_l=False)
        assert isinstance(k1, np.ndarray)
        expected_shape = (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
        assert (
            k1.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {k1.shape}"

    def test_knitting_matrix_computation_q0_raises_error(self, fibonacci_model: AnyonModel) -> None:
        """Test that knitting matrix computation for q=0 raises AssertionError as observed in notebook."""
        with pytest.raises(AssertionError, match="q must be strictly positive"):
            fibonacci_model.compute_knitting_matrix(q=0)

    def test_fusion_matrix_properties(self, fibonacci_model: AnyonModel) -> None:
        """Test properties of the fusion matrix."""
        fusion_matrix = fibonacci_model.n_symbols

        # Test shape
        assert fusion_matrix.shape == (2, 2, 2)

        # Test specific fusion rules as defined in notebook
        # Identity (0) x Identity (0) = Identity (0)
        assert fusion_matrix[0, 0, 0] == 1
        assert fusion_matrix[0, 0, 1] == 0

        # Identity (0) x Tau (1) = Tau (1)
        assert fusion_matrix[0, 1, 1] == 1
        assert fusion_matrix[0, 1, 0] == 0

        # Tau (1) x Identity (0) = Tau (1)
        assert fusion_matrix[1, 0, 1] == 1
        assert fusion_matrix[1, 0, 0] == 0

        # Tau (1) x Tau (1) can give both Identity (0) and Tau (1)
        assert fusion_matrix[1, 1, 0] == 1
        assert fusion_matrix[1, 1, 1] == 1

    def test_f_matrix_properties(self, fibonacci_model: AnyonModel) -> None:
        """Test properties of the F matrix."""
        f_matrix = fibonacci_model.f_matrix

        # Test shape
        assert f_matrix.shape == (2, 2, 2, 2, 2, 2)

        # Test that it's complex
        assert np.iscomplexobj(f_matrix)

    def test_r_matrix_properties(self, fibonacci_model: AnyonModel) -> None:
        """Test properties of the R matrix."""
        r_matrix = fibonacci_model.r_matrix

        # Test shape
        assert r_matrix.shape == (2, 2, 2)

        # Test that it's complex
        assert np.iscomplexobj(r_matrix)

    def test_braiding_matrix_properties(self, fibonacci_model: AnyonModel) -> None:
        """Test properties of the braiding matrix."""
        b_matrix = fibonacci_model.b_matrix

        # Test that it's complex
        assert np.iscomplexobj(b_matrix)

        # Test that it's not all zeros
        assert not np.allclose(b_matrix, 0)

    def test_check_rule_method(self, fibonacci_model: AnyonModel) -> None:
        """Test the check_rule method with various anyon combinations."""
        # Test valid fusion rules
        assert fibonacci_model.check_rule(np.array([0]), np.array([0]), np.array([0]))[
            0
        ]
        assert fibonacci_model.check_rule(np.array([0]), np.array([1]), np.array([1]))[
            0
        ]
        assert fibonacci_model.check_rule(np.array([1]), np.array([0]), np.array([1]))[
            0
        ]
        assert fibonacci_model.check_rule(np.array([1]), np.array([1]), np.array([0]))[
            0
        ]
        assert fibonacci_model.check_rule(np.array([1]), np.array([1]), np.array([1]))[
            0
        ]

        # Test invalid fusion rules
        assert not fibonacci_model.check_rule(
            np.array([0]), np.array([0]), np.array([1])
        )[0]
        assert not fibonacci_model.check_rule(
            np.array([0]), np.array([1]), np.array([0])
        )[0]
        assert not fibonacci_model.check_rule(
            np.array([1]), np.array([0]), np.array([0])
        )[0]

    def test_complex_calculations_dont_crash(self, fibonacci_model: AnyonModel) -> None:
        """Test that complex matrix calculations don't crash."""
        # These should all complete without error
        b_matrix = fibonacci_model.b_matrix
        l1 = fibonacci_model._compute_l_matrix(q=1)
        k1 = fibonacci_model.compute_knitting_matrix(q=1, return_l=False)
        assert isinstance(k1, np.ndarray)

        # Basic sanity checks
        assert not np.any(np.isnan(b_matrix))
        assert not np.any(np.isnan(l1))
        assert not np.any(np.isnan(k1))

        assert not np.any(np.isinf(b_matrix))
        assert not np.any(np.isinf(l1))
        assert not np.any(np.isinf(k1))


if __name__ == "__main__":
    pytest.main([__file__])
