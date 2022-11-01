import os

import numpy as np
import pytest

from tqsim import AnyonicCircuit, generate_braiding_operator, generate_basis


def test_init_1():
    circuit = AnyonicCircuit()

    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudits == 3
    assert circuit.dim == 3


def test_init_2():
    circuit = AnyonicCircuit(1, 4)

    assert circuit
    assert circuit.nb_qudits == 1
    assert circuit.nb_anyons_per_qudits == 4
    assert circuit.dim == 5


def test_init_3():
    circuit = AnyonicCircuit(2, 4)

    assert circuit
    assert circuit.nb_qudits == 2
    assert circuit.nb_anyons_per_qudits == 4
    assert circuit.dim == 34


def test_save():
    circuit = AnyonicCircuit()
    config_path = os.path.join(os.path.expanduser("~"), f".tqsim")
    store_path = os.path.join(config_path, "store")
    data_path = os.path.join(store_path, "1_3")
    basis_path = os.path.join(data_path, "basis.dat")
    sigmas_path = os.path.join(data_path, "sigmas.dat")

    assert os.path.exists(store_path)
    assert os.path.exists(data_path)
    assert os.path.exists(basis_path)
    assert os.path.exists(sigmas_path)


def test_initialize_1():
    circuit = AnyonicCircuit()
    try:
        circuit.initialize(np.ones(3) / np.sqrt(3))
    except Exception:
        assert False


def test_initialize_2():
    circuit = AnyonicCircuit()
    with pytest.raises(ValueError):
        circuit.initialize(np.ones(3))


def test_initialize_3():
    circuit = AnyonicCircuit()
    with pytest.raises(ValueError):
        circuit.initialize(np.ones(5))


def test_initialize_4():
    circuit = AnyonicCircuit()
    circuit.braid(1, 2)
    with pytest.raises(Exception):
        circuit.initialize(np.ones(3) / np.sqrt(3))


def test_braid_1():
    circuit = AnyonicCircuit()
    try:
        circuit.braid(1, 2)
        circuit.braid(2, 3)
    except:
        assert False


def test_braid_2():
    circuit = AnyonicCircuit()
    try:
        circuit.initialize(np.ones(3) / np.sqrt(3))
        circuit.braid(1, 2)
        circuit.braid(2, 3)
    except:
        assert False


def test_braid_3():
    circuit = AnyonicCircuit()
    circuit.measure()

    with pytest.raises(Exception):
        circuit.braid(1, 2)


def test_braid_4():
    circuit = AnyonicCircuit()

    with pytest.raises(Exception):
        circuit.braid(1, 3)


def test_braid_5():
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(0, 1)


def test_braid_6():
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(1, 0)


def test_braid_7():
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(4, 3)


def test_braid_8():
    circuit = AnyonicCircuit()

    with pytest.raises(ValueError):
        circuit.braid(3, 4)


def test_measure_1():
    circuit = AnyonicCircuit()
    try:
        circuit.measure()
    except:
        assert False
