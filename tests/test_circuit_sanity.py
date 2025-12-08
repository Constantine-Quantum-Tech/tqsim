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

from tqsim import AnyonicCircuit
from tqsim.tools.braiding_generators.fib_qudit import braiding_generator
from tqsim.tools.braiding_generators.fib_multi_qudits import (
    braiding_generator as multi_fib_braiding_generator,
)


def test_hadamard_gate_sanity() -> None:
    """Test Hadamard gate implementation.

    This test verifies:
    1. The Hadamard matrix from quantum circuit matches the one from braiding_generator
    2. The matrix is unitary (U @ U.conj().T = I)
    3. H^2 is close to identity (up to small error)
    """
    # Hadamard weaving sequence
    hadamard_seq = [
        [1, -4],
        [2, 2],
        [1, -2],
        [2, 2],
        [1, -2],
        [2, -2],
        [1, 2],
        [2, -4],
        [1, -2],
        [2, 2],
        [1, 2],
        [2, -2],
        [1, -2],
    ]

    # Create circuit and apply Hadamard sequence
    circuit = AnyonicCircuit(nb_qudits=1, nb_anyons_per_qudit=3)
    circuit.braid_sequence(hadamard_seq)

    # Get the unitary from the circuit
    hadamard_circuit_unitary = circuit.unitary()

    # Compute the Hadamard matrix using braiding_generator tool
    # For a single qudit with 3 anyons, we need to compute the product of braiding operations
    # The braiding sequence uses operators sigma_1 and sigma_2
    sigma_1, basis = braiding_generator(index=1, n_anyons=3, show=False)
    sigma_2, _ = braiding_generator(index=2, n_anyons=3, show=False)

    sigma_1 = np.array(sigma_1)
    sigma_2 = np.array(sigma_2)

    # Apply the braiding sequence using matrix multiplication
    hadamard_generator_unitary = np.eye(3, dtype=complex)
    for op_index, power in hadamard_seq:
        if op_index == 1:
            operator = sigma_1
        else:  # op_index == 2
            operator = sigma_2

        if power > 0:
            for _ in range(power):
                hadamard_generator_unitary = operator @ hadamard_generator_unitary
        else:  # negative power means inverse
            operator_inv = np.linalg.inv(operator)
            for _ in range(abs(power)):
                hadamard_generator_unitary = operator_inv @ hadamard_generator_unitary

    # Test 1: Compare circuit unitary with braiding_generator result
    # They should be equal up to numerical precision
    for s1, state_1 in enumerate(circuit.basis):
        # {'qudits': [[i, j], [k, l]], 'roots': [m]}
        tested_state_1 = [int(state_1.charges[0]), int(state_1.charges[1])]
        idx_1 = basis.index(tested_state_1)
        for s2, state_2 in enumerate(circuit.basis):
            tested_state_2 = [int(state_2.charges[0]), int(state_2.charges[1])]
            idx_2 = basis.index(tested_state_2)
            assert np.isclose(
                hadamard_circuit_unitary[s1, s2],
                hadamard_generator_unitary[idx_1, idx_2],
                rtol=1e-10,
            ), "Circuit Hadamard matrix does not match braiding_generator result"

    # Test 2: Check unitarity (U @ U.conj().T = I)
    identity_check = hadamard_circuit_unitary @ hadamard_circuit_unitary.conj().T
    assert np.allclose(
        identity_check, np.eye(3), atol=1e-10
    ), "Hadamard matrix is not unitary"


def test_cnot_gate_sanity() -> None:
    """Test Hadamard gate implementation.

    This test verifies:
    1. The Hadamard matrix from quantum circuit matches the one from braiding_generator
    2. The matrix is unitary (U @ U.conj().T = I)
    3. H^2 is close to identity (up to small error)
    """
    # CNOT weaving sequence
    cnot_seq = [
        [3, 1],
        [4, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [1, -1],
        [2, -1],
        [2, -1],
        [1, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [1, -1],
        [2, -1],
        [2, -1],
        [1, -1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [1, -1],
        [2, -1],
        [2, -1],
        [1, -1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [1, 1],
        [2, 1],
        [2, 1],
        [1, 1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [1, -1],
        [2, -1],
        [2, -1],
        [1, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [2, -1],
        [3, -1],
        [3, -1],
        [2, -1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [4, 1],
        [3, 1],
        [3, 1],
        [4, 1],
        [2, 1],
        [3, 1],
        [3, 1],
        [2, 1],
        [4, -1],
        [3, -1],
        [3, -1],
        [4, -1],
        [4, -1],
        [3, -1],
    ]

    # Create circuit and apply Hadamard sequence
    circuit = AnyonicCircuit(nb_qudits=2, nb_anyons_per_qudit=3)
    circuit.braid_sequence(cnot_seq)

    # Get the unitary from the circuit
    cnot_circuit_unitary = circuit.unitary()

    # Compute the CNOT matrix using braiding_generator tool
    sigmas = []
    for s in range(5):
        sigma, basis = multi_fib_braiding_generator(
            index=s + 1, n_qudits=2, qudit_len=2, show=False
        )
        sigmas.append(np.array(sigma, dtype=np.complex128))

    # Apply the braiding sequence using matrix multiplication
    cnot_generator_unitary = np.eye(len(basis), dtype=complex)
    for op_index, power in cnot_seq:
        operator = sigmas[op_index - 1]
        if power > 0:
            for _ in range(power):
                cnot_generator_unitary = operator @ cnot_generator_unitary
        else:  # negative power means inverse
            operator_inv = np.linalg.inv(operator)
            for _ in range(abs(power)):
                cnot_generator_unitary = operator_inv @ cnot_generator_unitary

    # Test 1: Compare circuit unitary with braiding_generator result
    # They should be equal up to numerical precision
    for s1, state_1 in enumerate(circuit.basis):
        # {'qudits': [[i, j], [k, l]], 'roots': [m]}
        tested_state_1 = {
            "qudits": [
                [int(state_1.charges[q * 2]), int(state_1.charges[q * 2 + 1])]
                for q in range(circuit.nb_qudits)
            ],
            "roots": [
                int(state_1.charges[a])
                for a in range(circuit.nb_qudits * 2, len(state_1.charges))
            ],
        }
        idx_1 = basis.index(tested_state_1)
        for s2, state_2 in enumerate(circuit.basis):
            tested_state_2 = {
                "qudits": [
                    [int(state_2.charges[q * 2]), int(state_2.charges[q * 2 + 1])]
                    for q in range(circuit.nb_qudits)
                ],
                "roots": [
                    int(state_2.charges[a])
                    for a in range(circuit.nb_qudits * 2, len(state_2.charges))
                ],
            }
            idx_2 = basis.index(tested_state_2)
            assert np.isclose(
                cnot_circuit_unitary[s1, s2],
                cnot_generator_unitary[idx_1, idx_2],
                rtol=1e-10,
            )
        "Circuit Hadamard matrix does not match braiding_generator result"

    # Test 2: Check unitarity (U @ U.conj().T = I)
    identity_check = cnot_circuit_unitary @ cnot_circuit_unitary.conj().T
    assert np.allclose(
        identity_check, np.eye(len(basis)), atol=1e-10
    ), "Hadamard matrix is not unitary"
