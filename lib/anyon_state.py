from abc import ABC, abstractmethod
import numpy as np


class AnyonState(ABC):
    """Abstract base class for anyonic states."""
    @abstractmethod
    def is_valid(self, model):
        pass

class StandardAnyonState(AnyonState):
    r"""Anyonic State involving multiple anyons in a linear 
    fusion tree (standard basis).
    
    Example:

    i0  i1  i2  i3
         \ /   /    /
          o1  /    /
           \ /    /
            o2   /
             \  /
              o3

        inputs:
        -------
            inputs   : np.array([i0, i1, i2, ..])
            outcomes : np.array([o1, o2, o3, ..])
    """
    def __init__(self, inputs: np.ndarray, outcomes: np.ndarray):
        self.inputs = inputs
        self.outcomes = outcomes
        assert len(outcomes) == len(inputs) - 1, "Number of outcomes must be one less than number of inputs."

    def is_valid(self, model):
        return model.check_state(self)

class SparseAnyonState(AnyonState):
    r"""Anyonic State involving multiple qudits, 
    each qudit being made of multiple anyons (sparse basis).
    
    Example:
        a b ..        g  h
        \/  / \/  / \/  /
        i\ /  k\ /  m\ /
          \     /     /
          j\  l/     /n
            \ /     /
            s\     /
              \   /
               \ /
               t|
        Here, 
            charges = [a, b, .., g, h, i, j, k, l, m, n, s, t]
            nb_qudits = 3
            nb_anyons_per_qudit = 3
    """
    def __init__(self, charges: np.ndarray, nb_qudits: int, nb_anyons_per_qudit: int):
        self.charges = charges
        self.nb_qudits = nb_qudits
        self.nb_anyons_per_qudit = nb_anyons_per_qudit
        assert len(charges) == nb_qudits * (2 * nb_anyons_per_qudit - 1) + (nb_qudits - 1), "Charges length does not match the number of qudits and anyons per qudit."

    def is_valid(self, model):
        return model.check_state(self)

class ComputationalSparseAnyonState(AnyonState):
    r"""Anyonic State involving multiple qudits, 
    each qudit being made of multiple anyons (sparse basis)
    such that the input anyons have identical charges.

    Example:
        a a ..        a  a
        \/  / \/  / \/  /
        i\ /  k\ /  m\ /
          \     /     /
          j\  l/     /n
            \ /     /
            s\     /
              \   /
               \ /
               t|
        Here, 
            charges = [i, j, k, l, m, n, s, t]
            nb_qudits = 3
            nb_anyons_per_qudit = 3
    """
    def __init__(self, charges: np.ndarray, nb_qudits: int, nb_anyons_per_qudit: int, input_charge: int):
        self.charges = charges
        self.nb_qudits = nb_qudits
        self.nb_anyons_per_qudit = nb_anyons_per_qudit
        self.input_charge = input_charge
        assert len(charges) == nb_qudits * (nb_anyons_per_qudit - 1) + (nb_qudits - 1), "Charges length does not match the number of qudits and anyons per qudit."

    def is_valid(self, model):
        sparse_state = SparseAnyonState(
            charges = np.concatenate((
                np.full(self.nb_qudits * self.nb_anyons_per_qudit, self.input_charge),
                self.charges
            )),
            nb_qudits = self.nb_qudits,
            nb_anyons_per_qudit = self.nb_anyons_per_qudit
        )
        return model.check_state(sparse_state)
