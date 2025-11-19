import itertools
from typing import Any

import numpy as np

from tqsim.models.fibonacci import FIBONACCI_MODEL
from tqsim.models.ising import ISING_MODEL
from tqsim.tools.braiding_generators.fib_multi_qudits import (
    braiding_matrix as fib_b_matrix,
    f_matrix as fib_f_matrix,
    knitting_matrix as fib_k_matrix,
    l_matrix as fib_l_matrix,
)
from tqsim.tools.braiding_generators.ising_multi_qudits import (
    braiding_matrix as ising_b_matrix,
    knitting_matrix as ising_k_matrix,
    l_matrix as ising_l_matrix,
)

fib_model = FIBONACCI_MODEL
ising_model = ISING_MODEL


def test_f_matrix() -> None:
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = fib_model.f_matrix[a, b, c, d, :, :]
        got = fib_f_matrix(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_b_matrix() -> None:
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = fib_model.b_matrix[a, b, c, d, :, :]
        got = fib_b_matrix(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_ising_b_matrix() -> None:
    # iterate over all possible a, b, c, d in [0, 1]
    for a, b, c, d in itertools.product([0, 1], repeat=4):
        expected = ising_model.b_matrix[a, b, c, d, :, :]
        got = ising_b_matrix(a, b, c, d)
        assert np.isclose(got, expected, atol=1e-6).all()


def test_l_matrix_fib_and_ising() -> None:
    # Check L matrix mapping for both Fibonacci and Ising models for q=3
    q = 3
    a = [1 for _ in range(q + 2)]

    # Fibonacci
    l_matrix = fib_model._compute_l_matrix(q)
    for params in itertools.product([0, 1], repeat=2 + 2 + 2 * q):
        k, h, i_, i = params[0], params[1], params[2], params[3]
        jj_ = list(params[4 : 4 + q])
        jj = list(params[4 + q : 4 + 2 * q])

        idx = tuple(a + [h, k, i] + jj + [i_] + jj_)
        l1 = l_matrix[idx]
        l2 = fib_l_matrix(k, h, i_, i, jj_, jj)
        assert np.isclose(l1, l2), f"Fibonacci L mismatch for params: {params}"

    # Ising
    l_matrix = ising_model._compute_l_matrix(q)
    for params in itertools.product([0, 1], repeat=2 + 2 + 2 * q):
        k, h, i_, i = params[0], params[1], params[2], params[3]
        jj_ = list(params[4 : 4 + q])
        jj = list(params[4 + q : 4 + 2 * q])

        idx = tuple(a + [h, k, i] + jj + [i_] + jj_)
        l1 = l_matrix[idx]
        l2 = ising_l_matrix(k, h, i_, i, jj_, jj)
        assert np.isclose(l1, l2), f"Ising L mismatch for params: {params}"


def test_k_matrix_fib_and_ising() -> None:
    # Check K (knitting) matrix mapping for both Fibonacci
    # and Ising models for q=3
    q = 2
    a = [1 for _ in range(q + 2)]

    # Fibonacci
    k_matrix = fib_model.compute_knitting_matrix(q)

    def new_k(
        jm: int, jmo: int, jmoo: int, jmo_: int, h: int, i_: int, i: int, jj_: list[int], jj: list[int]
    ) -> Any:
        return k_matrix[tuple(a + [h, jmoo, jm, jmo, i] + jj + [jmo_, i_] + jj_)]  # type: ignore[call-overload]

    for params in itertools.product([0, 1], repeat=3 + 3 + 2 + 2 * q):
        jm, jmo, jmoo, jmo_ = params[0], params[1], params[2], params[3]
        h, i_, i = params[4], params[5], params[6]
        jj_ = list(params[7 : 7 + q])
        jj = list(params[7 + q : 7 + 2 * q])

        k1 = new_k(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        k2 = fib_k_matrix(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        assert np.isclose(k1, k2), f"Fibonacci K mismatch for params: {params}"

    # Ising
    k_matrix = ising_model.compute_knitting_matrix(q)

    def new_k_ising(
        jm: int, jmo: int, jmoo: int, jmo_: int, h: int, i_: int, i: int, jj_: list[int], jj: list[int]
    ) -> Any:
        return k_matrix[tuple(a + [h, jmoo, jm, jmo, i] + jj + [jmo_, i_] + jj_)]  # type: ignore[call-overload]

    for params in itertools.product([0, 1], repeat=3 + 3 + 2 + 2 * q):
        jm, jmo, jmoo, jmo_ = params[0], params[1], params[2], params[3]
        h, i_, i = params[4], params[5], params[6]
        jj_ = list(params[7 : 7 + q])
        jj = list(params[7 + q : 7 + 2 * q])

        k1 = new_k_ising(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        k2 = ising_k_matrix(jm, jmo, jmoo, jmo_, h, i_, i, jj_, jj)
        assert np.isclose(k1, k2), f"Ising K mismatch for params: {params}"
