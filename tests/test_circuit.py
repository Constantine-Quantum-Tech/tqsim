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

import itertools
import os

import numpy as np
import pytest

from tqsim import AnyonicCircuit
from tqsim.lib.anyon_model import AnyonModel
from tqsim.models.fibonacci import FIBONACCI_MODEL
from tqsim.models.ising import ISING_MODEL


def clean_cache() -> None:
    config_path = os.path.join(os.path.expanduser("~"), ".tqsim")
    store_path = os.path.join(config_path, "store")
    temp_path = os.path.join(config_path, "temp_store")

    if os.path.exists(store_path):
        os.rename(store_path, temp_path)


def pytest_sessionstart(session: pytest.Session) -> None:
    """Run before any tests are collected or executed."""
    print("Running pre-test setup...")
    clean_cache()


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    """Called after the whole test run completes."""
    print("\nAll tests finished. Cleaning up...")
    # remove temporary stored files after tests
    config_path = os.path.join(os.path.expanduser("~"), ".tqsim")
    store_path = os.path.join(config_path, "store")
    temp_path = os.path.join(config_path, "temp_store")

    if os.path.exists(store_path):
        os.rmdir(store_path)
    if os.path.exists(temp_path):
        os.rename(temp_path, store_path)


def test_init_1() -> None:
    circuit = AnyonicCircuit()

    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudit == 3
    assert circuit.dim == 3
    assert circuit.braiding_operators
    assert circuit.basis
    assert len(circuit.braiding_operators) == 2
    assert circuit.braiding_operators[0].shape == (3, 3)


def test_init_2() -> None:
    circuit = AnyonicCircuit(1, 4)

    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudit == 4
    assert circuit.dim == 5


def test_init_3() -> None:
    circuit = AnyonicCircuit(2, 4)

    assert circuit
    assert circuit.nb_qudits == 2
    assert circuit.nb_anyons_per_qudit == 4
    assert circuit.dim == 34


def test_init_model_1() -> None:
    # Z_N model (Abelian model)
    num_charges = 5
    n_symbols = np.zeros((num_charges, num_charges, num_charges), dtype=int)
    for i, j, k in itertools.product(range(num_charges), repeat=3):
        if (i + j) % num_charges == k:
            n_symbols[i, j, k] = 1

    f_matrix = np.zeros(
        (num_charges, num_charges, num_charges, num_charges, num_charges, num_charges),
        dtype=complex,
    )
    for i, j, k, total_charge, m, n in itertools.product(range(num_charges), repeat=6):
        if (
            (i + j + k) % num_charges == total_charge
            and (i + j) % num_charges == m
            and (j + k) % num_charges == n
        ):
            f_matrix[i, j, k, total_charge, m, n] = 1

    r_matrix = np.zeros((num_charges, num_charges, num_charges), dtype=complex)
    for i, j, k in itertools.product(range(num_charges), repeat=3):
        if (i + j) % num_charges == k:
            r_matrix[i, j, k] = np.exp(2j * np.pi * i * j / num_charges)

    zn_model = AnyonModel(n_symbols, f_matrix, r_matrix, name="Z_N")

    circuit = AnyonicCircuit(
        nb_qudits=1, nb_anyons_per_qudit=4, model=zn_model, input_charge=1
    )

    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudit == 4
    assert circuit.input_charge == 1  # Which charge is used for braiding
    assert circuit.model.name == "Z_N"
    assert circuit.dim == 1
    assert len(circuit.braiding_operators) == 3
    assert circuit.braiding_operators[0].shape == (1, 1)

    try:
        circuit.model = FIBONACCI_MODEL  # type: ignore[misc]
        assert False
    except AttributeError:
        assert True


def test_init_model_2() -> None:
    circuit = AnyonicCircuit(
        nb_qudits=1, nb_anyons_per_qudit=3, model=ISING_MODEL, input_charge=1
    )
    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudit == 3
    assert circuit.input_charge == 1
    assert circuit.model.name == "Ising"
    assert circuit.dim == 2
    assert len(circuit.braiding_operators) == 2
    assert circuit.braiding_operators[0].shape == (2, 2)

    try:
        circuit.model = FIBONACCI_MODEL  # type: ignore[misc]
        assert False
    except AttributeError:
        assert True


def test_save() -> None:
    _ = AnyonicCircuit()
    config_path = os.path.join(os.path.expanduser("~"), ".tqsim")
    store_path = os.path.join(config_path, "store")
    data_path = os.path.join(store_path, "Fibonacci-1-3-1")
    basis_path = os.path.join(data_path, "-basis.dat")
    sigmas_path = os.path.join(data_path, "-sigmas.dat")

    assert os.path.exists(store_path)
    assert os.path.exists(data_path)
    assert os.path.exists(basis_path)
    assert os.path.exists(sigmas_path)


def test_initialize_1() -> None:
    circuit = AnyonicCircuit()
    try:
        circuit.initialize(np.ones(3) / np.sqrt(3))
    except Exception:
        assert False


def test_initialize_2() -> None:
    circuit = AnyonicCircuit()
    with pytest.raises(ValueError):
        circuit.initialize(np.ones(3))


def test_initialize_3() -> None:
    circuit = AnyonicCircuit()
    with pytest.raises(ValueError):
        circuit.initialize(np.ones(5))


def test_initialize_4() -> None:
    circuit = AnyonicCircuit()
    circuit.braid(1, 2)
    with pytest.raises(Exception):
        circuit.initialize(np.ones(3) / np.sqrt(3))


def test_braid_1() -> None:
    circuit = AnyonicCircuit()
    try:
        circuit.braid(1, 2)
        circuit.braid(2, 3)
    except Exception:
        assert False


def test_braid_2() -> None:
    circuit = AnyonicCircuit()
    try:
        circuit.initialize(np.ones(3) / np.sqrt(3))
        circuit.braid(1, 2)
        circuit.braid(2, 3)
    except Exception:
        assert False


def test_braid_3() -> None:
    circuit = AnyonicCircuit()
    circuit.measure()

    with pytest.raises(Exception):
        circuit.braid(1, 2)


def test_braid_4() -> None:
    circuit = AnyonicCircuit()

    with pytest.raises(Exception):
        circuit.braid(1, 3)


def test_braid_5() -> None:
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(0, 1)


def test_braid_6() -> None:
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(1, 0)


def test_braid_7() -> None:
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(4, 3)


def test_braid_8() -> None:
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(3, 4)


def test_measure_1() -> None:
    circuit = AnyonicCircuit()
    try:
        circuit.measure()
    except Exception:
        assert False
