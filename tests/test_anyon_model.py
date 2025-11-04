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

from tqsim.lib.anyon_model import AnyonModel


class TestAnyonModel:
    """Test suite for AnyonModel class based on notebook experiments."""

    @pytest.fixture
    def fibonacci_model(self):
        """Create a Fibonacci anyon model as tested in the notebook."""
        # Fusion matrix setup from notebook
        fusion_matrix = np.zeros((2, 2, 2))
        for a1 in range(2):
            for a2 in range(2):
                for outcome in range(2):
                    if (a1, a2) == (1, 1):
                        fusion_matrix[a1, a2, outcome] = 1
                    else:
                        if (a1 + a2) == outcome:
                            fusion_matrix[a1, a2, outcome] = 1

        def F(a1, a2, a3, outcome):
            """F matrix helper function from notebook."""
            inv_phi = (np.sqrt(5) - 1) / 2  # inverse of golden number
            f_matrix = np.array([[0, 0], [0, 0]])

            # a1 + a2 + a3 + outcome = 4
            if a1 + a2 + a3 + outcome == 4:
                f_matrix = np.array(
                    [[inv_phi, np.sqrt(inv_phi)], [np.sqrt(inv_phi), -inv_phi]]
                )
            # a1 + a2 + a3 + outcome = 3
            elif a1 + a2 + a3 + outcome == 3:
                f_matrix = np.array([[0, 0], [0, 1]])
            # a1 + a2 + a3 + outcome = 2
            elif a1 + a2 + a3 + outcome == 2:
                if a1 + a2 == 2:
                    f_matrix = np.array([[0, 1], [0, 0]])
                elif a2 + a3 == 2:
                    f_matrix = np.array([[0, 0], [1, 0]])
                elif a1 + a3 == 2:
                    f_matrix = np.array([[0, 0], [0, 1]])
                elif a3 + outcome == 2:
                    f_matrix = np.array([[0, 1], [0, 0]])
                elif a1 + outcome == 2:
                    f_matrix = np.array([[0, 0], [1, 0]])
                elif a2 + outcome == 2:
                    f_matrix = np.array([[0, 0], [0, 1]])
            # a1 + a2 + a3 + outcome = 0
            elif a1 + a2 + a3 + outcome == 0:
                f_matrix = np.array([[1, 0], [0, 0]])

            return f_matrix

        def R(a1, a2):
            """R matrix helper function from notebook."""
            if a1 + a2 == 2:
                r_matrix = np.array(
                    [
                        [np.exp(-4 * np.pi * 1j / 5), 0],
                        [0, np.exp(3 * np.pi * 1j / 5)],
                    ]
                )
            else:
                r_matrix = np.array([[1, 0], [0, 1]])
            return r_matrix

        # Build F and R matrices
        F_matrix = np.zeros((2, 2, 2, 2, 2, 2)) * (1 + 0j)
        R_matrix = np.zeros((2, 2, 2)) * (1 + 0j)

        for a1 in range(2):
            for a2 in range(2):
                for a3 in range(2):
                    for outcome in range(2):
                        F_matrix[a1, a2, a3, outcome] = F(a1, a2, a3, outcome)

        for a1 in range(2):
            for a2 in range(2):
                R_matrix[a1, a2] = R(a1, a2).diagonal()

        return AnyonModel(fusion_matrix, F_matrix, R_matrix)

    def test_model_initialization(self, fibonacci_model):
        """Test that the model initializes correctly."""
        assert fibonacci_model is not None
        assert hasattr(fibonacci_model, "N_symbols")
        assert hasattr(fibonacci_model, "F_matrix")
        assert hasattr(fibonacci_model, "R_matrix")
        assert hasattr(fibonacci_model, "B_matrix")

    def test_braiding_matrix_shape(self, fibonacci_model):
        """Test that braiding matrix has correct shape as observed in notebook."""
        B = fibonacci_model.B_matrix
        expected_shape = (2, 2, 2, 2, 2, 2)
        assert (
            B.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {B.shape}"

    def test_L_matrix_computation_q1(self, fibonacci_model):
        """Test L matrix computation for q=1 (successful case from notebook)."""
        L1 = fibonacci_model._compute_L_matrix(q=1)
        expected_shape = (2, 2, 2, 2, 2, 2, 2, 2, 2)
        assert (
            L1.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {L1.shape}"

    def test_L_matrix_computation_q0_raises_error(self, fibonacci_model):
        """Test that L matrix computation for q=0 raises AssertionError as observed in notebook."""
        with pytest.raises(AssertionError, match="q must be strictly positive"):
            fibonacci_model._compute_L_matrix(q=0)

    def test_knitting_matrix_computation_q1(self, fibonacci_model):
        """Test knitting matrix computation for q=1 (successful case from notebook)."""
        K1 = fibonacci_model.compute_knitting_matrix(q=1)
        expected_shape = (2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
        assert (
            K1.shape == expected_shape
        ), f"Expected shape {expected_shape}, got {K1.shape}"

    def test_knitting_matrix_computation_q0_raises_error(self, fibonacci_model):
        """Test that knitting matrix computation for q=0 raises AssertionError as observed in notebook."""
        with pytest.raises(AssertionError, match="q must be strictly positive"):
            fibonacci_model.compute_knitting_matrix(q=0)

    def test_fusion_matrix_properties(self, fibonacci_model):
        """Test properties of the fusion matrix."""
        fusion_matrix = fibonacci_model.N_symbols

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

    def test_F_matrix_properties(self, fibonacci_model):
        """Test properties of the F matrix."""
        F_matrix = fibonacci_model.F_matrix

        # Test shape
        assert F_matrix.shape == (2, 2, 2, 2, 2, 2)

        # Test that it's complex
        assert np.iscomplexobj(F_matrix)

    def test_R_matrix_properties(self, fibonacci_model):
        """Test properties of the R matrix."""
        R_matrix = fibonacci_model.R_matrix

        # Test shape
        assert R_matrix.shape == (2, 2, 2)

        # Test that it's complex
        assert np.iscomplexobj(R_matrix)

    def test_braiding_matrix_properties(self, fibonacci_model):
        """Test properties of the braiding matrix."""
        B = fibonacci_model.B_matrix

        # Test that it's complex
        assert np.iscomplexobj(B)

        # Test that it's not all zeros
        assert not np.allclose(B, 0)

    def test_check_rule_method(self, fibonacci_model):
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

    def test_complex_calculations_dont_crash(self, fibonacci_model):
        """Test that complex matrix calculations don't crash."""
        # These should all complete without error
        B = fibonacci_model.B_matrix
        L1 = fibonacci_model._compute_L_matrix(q=1)
        K1 = fibonacci_model.compute_knitting_matrix(q=1)

        # Basic sanity checks
        assert not np.any(np.isnan(B))
        assert not np.any(np.isnan(L1))
        assert not np.any(np.isnan(K1))

        assert not np.any(np.isinf(B))
        assert not np.any(np.isinf(L1))
        assert not np.any(np.isinf(K1))


if __name__ == "__main__":
    pytest.main([__file__])
