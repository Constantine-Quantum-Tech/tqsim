import unittest 
import numpy as np
from math import isclose
from cqt_anyons import gen_basis, braiding_generator

def spectral_distance(u : "np.array", w : "np.array"):
    """
    Error distance metric between SU(n) matrices using the Spectral Norm.
    returns 
    """
    u = u * (1 + 0j)
    w = w * (1 + 0j)
    global_phase_w = pow(np.cfloat(np.linalg.det(w)), 1/w.shape[0])
    global_phase_u = pow(np.cfloat(np.linalg.det(u)), 1/u.shape[0])

    diff = np.matrix(w / global_phase_w - u / global_phase_u)
    distance = abs((max(np.linalg.eig(diff.getH() @ diff)[0]))**(1/2))
    if distance > 1:
        return 2 - distance
    return distance

class FibonacciSequenceTester(unittest.TestCase): 
    
    def setUp(self): 
        pass
   
    def test_fibonacci_sequence(self):
        """
        Verifies wether the dimension of the fusion space obeys the Fibonacci
        series or not.

        Returns
        -------
        None.

        """
        fibo = [2, 3]
        for nb in range(1, 10):
            basis_size = len(gen_basis(nb_qudits=1, nb_anyons_per_qudit=2*nb))
            assert basis_size == fibo[0], "Fusion space's dimension of Fibonaci anyons should obey Fibonacci series."
            basis_size = len(gen_basis(nb_qudits=1, nb_anyons_per_qudit=2*nb+1))
            assert basis_size == fibo[1], "Fusion space's dimension of Fibonaci anyons should obey Fibonacci series."

            fibo[0] = sum(fibo)
            fibo[1] = sum(fibo)
        

class BraidingAlgebraTester(unittest.TestCase): 
    
    def setUp(self): 
        pass
   
    def test_yang_baxter_relations(self):
        """
        Verifies wether the braiding generators obey the Yang-Baxter equations
        or not.

        Returns
        -------
        None.

        """
        for nq in range(1, 3):
            for na in [3, 4]:
                nb_anyons = na * nq
                generators = [] 
                for index in range(1, nb_anyons):
                    generator = np.array(braiding_generator(index, nq, na))
                    generators.append(generator)
                    
                    # test unitarity
                    # si si^dagger = I
                    left = generator @ generator.T.conjugate()
                    identity = np.eye(generator.shape[0])
                    distance = spectral_distance(left, identity)
                    assert isclose(distance, 0, abs_tol=1e-15), f'sigma({index}, {nq}, {na}) is not unitary.'
                    
                for index_1 in range(0, nb_anyons - 1):
                    for index_2 in range(0, nb_anyons - 1):
                        if abs(index_1 - index_2) == 1:
                            # test second relation
                            # si sj si = sj si sj  /  |i-j| > 1
                            left = generators[index_1] @ generators[index_2] @ generators[index_1]
                            right = generators[index_2] @ generators[index_1] @ generators[index_2]
                            distance = spectral_distance(left, right)
                            assert isclose(distance, 0, abs_tol=1e-15), f'Yang Baxter Relation : {index_1} {nq} {na}.'
                        else:
                            # test third relation
                            # si sj = sj si  / |i-j| = 1
                            left = generators[index_1] @ generators[index_2]
                            right = generators[index_2] @ generators[index_1]
                            distance = spectral_distance(left, right)
                            assert isclose(distance, 0, abs_tol=1e-15), f'Commutative Yang Baxter relation : {index_1} {index_2} {nq} {na}.'

class NumericalTester(unittest.TestCase): 
    
    def setUp(self): 
        pass
    
    def test_hadamard(self):
        factor = 1 / np.sqrt(2)
        HADAMARD = np.array([[factor, factor], [factor, -factor]]) * 1j
        hadamard_sequence =[[1, -4], [2, 2], [1, -2], [2, 2], [1, -2], [2, -2],
                            [1, 2], [2, -4], [1, -2], [2, 2], [1, 2], [2, -2],
                            [1, -2]]
        
        for nq in range(1):
            for na in [3]:
                nb_anyons = nq * na
                generators = []
                for index in range(1, nb_anyons):
                    generator = np.array(braiding_generator(index, nq, na))
                    generators.append(generator)
                
                for qubit in range(nq):
                    m = np.eye(generator.shape[0])
                    for sigma, power in hadamard_sequence:
                        m = np.linalg.matrix_power(generators[sigma + qubit * na - 1], power) @ m
                    
                    basis = gen_basis(nq, na)
                    computational_m = np.zeros([2, 2]) * 1j
                    
                    row = 0
                    for f, basis_f in enumerate(basis):
                        if all(basis_f['qudits'][q][1] == 1 for q in range(nq)):
                            column = 0
                            for i, basis_i in enumerate(basis):
                                if all(basis_i['qudits'][q][1] == 1 for q in range(nq)):
                                    print(basis_f, basis_i)
                                    computational_m[row][column] = m[f][i]
                                    column += 1
                            
                            row += 1
                                
                    distance = spectral_distance(computational_m, HADAMARD)
                    print(distance)
                    print(computational_m)
                    assert isclose(distance, 0, abs_tol=0.00657)
    
if __name__ == '__main__': 
    unittest.main() 