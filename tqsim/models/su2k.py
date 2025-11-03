import itertools
import numpy as np

from tqsim.lib.anyon_model import AnyonModel


def get_fusion_matrix(k: int) -> np.ndarray:
    """
    Generates the SU(2)_k fusion matrix.

    Parameters
    ----------
    k : int
        The level of the SU(2)_k model.

    Returns
    -------
    np.ndarray
        A 3D numpy array representing the fusion matrix.
    """
    dim = k + 1
    fusion_matrix = np.zeros((dim, dim, dim), dtype=int)

    for a1 in range(dim):
        for a2 in range(dim):
            for outcome in range(dim):
                j1 = a1 / 2
                j2 = a2 / 2
                j = outcome / 2
                if abs(j1 - j2) <= j <= min(j1 + j2, k - j1 - j2):
                    fusion_matrix[a1, a2, outcome] = 1

    return fusion_matrix


def get_F_matrix(k: int):
    """
    Placeholder for the F matrix function for SU(2)_k.
    The actual implementation would depend on the specific details of the SU(2)_k model.

    Parameters
    ----------
    k : int
        The level of the SU(2)_k model.

    Returns
    -------
    function
        A function that computes the F matrix for given inputs.
    """
    q = np.exp(2j * np.pi / (k + 2))
    
    def crochet(n):
        """ 
        crochet = [n]_q 
                = (q**(n/2) - q**(-n/2)) / (q**(1/2) - q**(-1/2))
        """
        return (q**(n / 2) - q**(-n / 2)) / (q**(1 / 2) - q**(-1 / 2))

    def crochet_factorial(n):
        """ 
        [n]_q! = [n]_q [n-1]_q ... [1]_q
        """
        val = np.prod(crochet(np.arange(1, n + 1))) if n > 0 else 1
        # print(f"crochet_factorial({n}) = {val}")
        return val

    # delta = lambda j1, j2, j3: np.sqrt(
    def delta(j1, j2, j3):
        """ 
        Delta(j1, j2, j3) = sqrt(
            [-j1 + j2 + j3]_q! *
            [j1 - j2 + j3]_q! *
            [j1 + j2 - j3]_q! /
            [j1 + j2 + j3 + 1]_q!
        )
        """
        val =  np.sqrt(
        crochet_factorial(-j1 + j2 + j3) *
        crochet_factorial(j1 - j2 + j3) *
        crochet_factorial(j1 + j2 - j3) /
        crochet_factorial(j1 + j2 + j3 + 1)
    )
        return val

    def braceq(j1, j2, j3, j, j12, j23):
        """ 
        { j1 j2 j12 }
        { j3  j  j23 } = 
        delta(j1, j2, j12) *
        delta(j12, j3, j) *
        delta(j2, j3, j23) *
        delta(j1, j23, j) *
        sum over z of:
            (-1)^z *
            [z + 1]_q! /
            ( [z - j1 - j2 - j12]_q! *
              [z - j12 - j3 - j]_q! *
              [z - j2 - j3 - j23]_q! *
              [z - j1 - j23 - j]_q! *
              [j1 + j2 + j3 + j - z]_q! *
              [j1 + j12 + j3 + j23 - z]_q! *
              [j2 + j12 + j + j23 - z]_q! )
        """
        deltas = (
            delta(j1, j2, j12) * 
            delta(j12, j3, j) *
            delta(j2, j3, j23) *
            delta(j1, j23, j)
        )
        vals = []
        zmmin = int(max(j1 + j2 + j12, j12 + j3 + j, j2 + j3 + j23, j1 + j23 + j))
        zmmax = int(min(j1 + j2 + j3 + j, j1 + j12 + j3 + j23, j2 + j12 + j + j23))
        for z in range(zmmin, zmmax + 1):
            exp = (-1) ** z
            up =crochet_factorial(z + 1)
            down = (
                      crochet_factorial(z - j1 - j2 - j12)
                    * crochet_factorial(z - j12 - j3 - j)
                    * crochet_factorial(z - j2 - j3 - j)
                    * crochet_factorial(z - j2 - j3 - j23)
                    * crochet_factorial(j1 + j2 + j3 + j - z)
                    * crochet_factorial(j1 + j12 + j3 + j23 - z)
                    * crochet_factorial(j2 + j12 + j + j23 - z)
                )
            val = exp * up / down
            vals.append(val)

        sum_ = np.sum(vals)
        val = deltas * sum_
        return val

    F_matrix = np.zeros((k + 1, k + 1, k + 1, k + 1, k + 1, k + 1), dtype=complex)

    for a, b, c, d, e, f in itertools.product(range(k + 1), repeat=6):
        j1 = a / 2
        j2 = b / 2
        j3 = c / 2
        j  = d / 2
        j12= e / 2
        j23= f / 2

        F_matrix[a, b, c, d, e, f] = (
            (-1 + 0j) ** (j1 + j2 + j3 + j) * 
            braceq(j1, j2, j3, j, j12, j23) *
            np.sqrt(
                crochet(2 * j12 + 1) * crochet(2 * j23 + 1)
            )
        )

    return F_matrix


def get_R_matrix(k: int) -> np.ndarray:
    """
    Generates the SU(2)_k R matrix.

    Parameters
    ----------
    k : int
        The level of the SU(2)_k model.

    Returns
    -------
    np.ndarray
        A 3D numpy array representing the R matrix.
    """
    dim = k + 1
    R_matrix = np.zeros((dim, dim, dim), dtype=complex)
    q = np.exp(2j * np.pi / (k + 2))

    def crochet(n):
        """ 
        crochet = [n]_q 
                = (q**(n/2) - q**(-n/2)) / (q**(1/2) - q**(-1/2))
        """
        return (q**(n / 2) - q**(-n / 2)) / (q**(1 / 2) - q**(-1 / 2))
    
    for a, b, c in itertools.product(range(dim), repeat=3):
        j1 = a / 2
        j2 = b / 2
        j  = c / 2
        phase = q ** (
            (0.5 + 0j) * (j * (j + 1) - j1 * (j1 + 1) - j2 * (j2 + 1))
        )
        N_symbols = get_fusion_matrix(k)
        if N_symbols[a, b, c] != 0:
            R_matrix[a, b, c] = ((-1 + 0j) ** (j - j1 - j2)) * phase

    return R_matrix