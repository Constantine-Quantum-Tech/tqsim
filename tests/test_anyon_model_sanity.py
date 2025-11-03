import itertools
import numpy as np

from tools.braiding_generators.fib_multi_qudits import sigma as sigma_fib
from tools.braiding_generators.fib_multi_qudits import L as L_fib, B as B_fib
from tools.braiding_generators.fib_qudit import F as F_fib, B as B_fib_qudit
from tools.braiding_generators.fib_multi_qudits import S as K_fib

from tools.braiding_generators.ising_multi_qudits import sigma as sigma_ising
from tools.braiding_generators.ising_multi_qudits import (
    L as L_ising,
    B as B_ising,
)
from tools.braiding_generators.ising_qudit import (
    F as F_ising,
    B as B_ising_qudit,
)
from tools.braiding_generators.ising_multi_qudits import S as K_ising

from tqsim.models.fibonacci import FIBONACCI_MODEL
from tqsim.models.ising import ISING_MODEL


fib_model = FIBONACCI_MODEL
ising_model = ISING_MODEL


def test_F_matrix():
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = fib_model.F_matrix[a, b, c, d, :, :]
        got = F_fib(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_B_matrix():
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = fib_model.B_matrix[a, b, c, d, :, :]
        got = B_fib(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_ising_B_matrix():
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = ising_model.B_matrix[a, b, c, d, :, :]
        got = B_ising(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_L_matrix_fib_and_ising():
    # Check L matrix mapping for both Fibonacci and Ising models for q=3
    q = 3
    a = [1 for _ in range(q + 2)]

    # Fibonacci
    L_matrix = fib_model._compute_L_matrix(q)
    for params in itertools.product([0, 1], repeat=2 + 2 + 2 * q):
        k, h, i_, i = params[0], params[1], params[2], params[3]
        jj_ = list(params[4 : 4 + q])
        jj = list(params[4 + q : 4 + 2 * q])

        idx = tuple(a + [h, k, i] + jj + [i_] + jj_)
        l1 = L_matrix[idx]
        l2 = L_fib(k, h, i_, i, jj_, jj)
        assert np.isclose(l1, l2), f"Fibonacci L mismatch for params: {params}"

    # Ising
    L_matrix = ising_model._compute_L_matrix(q)
    for params in itertools.product([0, 1], repeat=2 + 2 + 2 * q):
        k, h, i_, i = params[0], params[1], params[2], params[3]
        jj_ = list(params[4 : 4 + q])
        jj = list(params[4 + q : 4 + 2 * q])

        idx = tuple(a + [h, k, i] + jj + [i_] + jj_)
        l1 = L_matrix[idx]
        l2 = L_ising(k, h, i_, i, jj_, jj)
        assert np.isclose(l1, l2), f"Ising L mismatch for params: {params}"


def test_K_matrix_fib_and_ising():
    # Check K (knitting) matrix mapping for both Fibonacci
    # and Ising models for q=3
    q = 2
    a = [1 for _ in range(q + 2)]

    # Fibonacci
    K_matrix = fib_model.compute_knitting_matrix(q)
    new_K = lambda jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj: K_matrix[
        tuple(a + [h, jmoo, jm, jmo, i] + jj + [jmo_, i_] + jj_)
    ]
    for params in itertools.product([0, 1], repeat=3 + 3 + 2 + 2 * q):
        jm, jmo, jmoo, jmo_ = params[0], params[1], params[2], params[3]
        h, i_, i = params[4], params[5], params[6]
        jj_ = list(params[7 : 7 + q])
        jj = list(params[7 + q : 7 + 2 * q])

        k1 = new_K(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        k2 = K_fib(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        assert np.isclose(k1, k2), f"Fibonacci K mismatch for params: {params}"

    # Ising
    K_matrix = ising_model.compute_knitting_matrix(q)
    new_K = lambda jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj: K_matrix[
        tuple(a + [h, jmoo, jm, jmo, i] + jj + [jmo_, i_] + jj_)
    ]
    for params in itertools.product([0, 1], repeat=3 + 3 + 2 + 2 * q):
        jm, jmo, jmoo, jmo_ = params[0], params[1], params[2], params[3]
        h, i_, i = params[4], params[5], params[6]
        jj_ = list(params[7 : 7 + q])
        jj = list(params[7 + q : 7 + 2 * q])

        k1 = new_K(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        k2 = K_ising(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        assert np.isclose(k1, k2), f"Ising K mismatch for params: {params}"
