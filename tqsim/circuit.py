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

import os
import pickle
from typing import List, Sequence, Tuple

import numpy as np

from .config import STORE_PATH  # For caching the bases and sigmas.
from .lib.basis_generator import generate_basis
from .lib.drawer import Drawer
from .lib.operator_generator import generate_braiding_operator


class AnyonicCircuit:
    """This class represents an anyon-based topological quantum circuit.
    Such a circuit is described by the number of qudits it contains, and
    the number of anyons each qudit contains.

    Example: >>>circuit = AnyonicCircuit(nb_qudits=2, nb_anyons_per_qudit=3)
                This circuit has 2 qudits, with 3 anyons each. The circuit has
                a total of 2*3=6 anyons.

    Parameters
    ----------
    nb_qudits : int, optional
        Number of qudits in the circuit. The default is 1.
    nb_anyons_per_qudit : int, optional
        Number of anyons in each qudit. The default is 3.

    Attributes
    ----------
    nb_qudits : int
        Number of qudits in the circuit.
    nb_anyons_per_qudit : int
        Number of anyons in each qudit.
    drawer : Drawer
        A drawer object that handles the drawing of the quantum circuit.
    dim : int
        The dimension of the fusion space for the quantum circuit.
    basis : List
        List of all the basis states.

    """

    def __init__(self, nb_qudits: int = 1, nb_anyons_per_qudit: int = 3):
        """
        Parameters
        ----------
        nb_qudits : int, optional
            Number of qudits in the circuit. The default is 1.
        nb_anyons_per_qudit : int, optional
            Number of anyons in each qudit. The default is 3.

        Returns
        -------
        None.

        """
        self.__nb_qudits = nb_qudits
        self.__nb_anyons_per_qudit = nb_anyons_per_qudit
        self.__nb_anyons = nb_qudits * nb_anyons_per_qudit

        self.__nb_braids: int = 0
        self.__braids_history: List[Tuple[int, int]] = []
        self.__measured: bool = False

        self.__basis, self.__dim = self.__get_basis()

        input_state = np.zeros((self.__dim, 1), dtype=np.complex128)
        input_state[0, 0] = 1
        self.__initial_state = input_state

        self.__sigmas = self.__get_sigmas()
        self.__unitary = np.eye(self.__dim)

        self.__drawer = Drawer(nb_qudits, nb_anyons_per_qudit)

    @property
    def nb_qudits(self):
        """Returns the number of qudits in the circuit.

        Returns
        -------
        int
            Number of qudits.
        """
        return self.__nb_qudits

    @property
    def nb_anyons_per_qudits(self):
        """Returns the number of anyons for each qudit in the circuit.

        Returns
        -------
        int
            Number of anyons per qudit.
        """
        return self.__nb_anyons_per_qudit

    @property
    def drawer(self):
        """Returns the drawer object for the circuit.

        Returns
        -------
        Drawer
            The circuit's drawer object.
        """
        return self.__drawer

    @property
    def dim(self):
        """Returns the dimension of the fusion space.

        Returns
        -------
        int
            Dimension of the fusion space.
        """
        return self.__dim

    @property
    def basis(self):
        """Returns a list of all the basis states for the circuit.

        Returns
        -------
        List
            List of all the basis states.

        """
        return self.__basis

    @property
    def braiding_operators(self):
        """Returns a list of all the braiding operators.

        Returns
        -------
        List
            List of all the braiding operators.

        """
        return self.__sigmas

    def __get_basis(self) -> Tuple[np.ndarray, int]:
        folder_path = os.path.join(
            STORE_PATH, f"{self.__nb_qudits}_{self.__nb_anyons_per_qudit}"
        )
        filename = os.path.join(folder_path, "basis.dat")
        try:
            with open(filename, "rb") as f:
                basis = pickle.load(f)
        except FileNotFoundError:
            basis = generate_basis(self.__nb_qudits, self.__nb_anyons_per_qudit)
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            with open(filename, "wb") as f:
                pickle.dump(basis, f)

        return basis, len(basis)

    def __get_sigmas(self) -> List[np.ndarray]:

        folder_path = os.path.join(
            STORE_PATH, f"{self.__nb_qudits}_{self.__nb_anyons_per_qudit}"
        )
        filename = os.path.join(folder_path, "sigmas.dat")
        try:
            with open(filename, "rb") as f:
                sigmas = pickle.load(f)
        except FileNotFoundError:
            sigmas = []
            for index in range(1, self.__nb_anyons):
                sigma = generate_braiding_operator(
                    index, self.__nb_qudits, self.__nb_anyons_per_qudit
                )
                sigmas.append(np.array(sigma))

            os.makedirs(os.path.dirname(filename), exist_ok=True)
            with open(filename, "wb") as f:
                pickle.dump(sigmas, f)

        return sigmas

    def initialize(self, input_state: np.ndarray):
        """Initializes the circuit in the state input_state.

        Parameters
        ----------
        input_state : np.ndarray
            A normalized quantum state with the same dimensions as the
            fusion space.
            Example:    for a 1-qudit circuit with 3 anyons, input_state must
                        be a 3 dimensional vector with norm 1.

        Raises
        ------
        Exception
            Will be raised if an initialization is attempted after performing
            braiding operations.
        ValueError
            Will be raised if the input state has the wrong dimension or
            is not normalized.

        Returns
        -------
        AnyonicCircuit
            A reference to the same circuit.
        """
        if self.__nb_braids > 0:
            raise Exception(
                "Initialization should happen before any braiding operation is performed!"
            )

        input_state = np.array(input_state)
        if not np.size(input_state) == self.__dim:
            raise ValueError(f"The state has wrong dimension. Should be {self.__dim}")

        norm = np.sum(np.real(input_state * input_state.conjugate()))
        if not np.isclose(norm, 1):
            raise ValueError("The input state is not normalized correctly!")

        self.__initial_state = np.reshape(input_state, (self.__dim, 1))

        return self

    def braid(self, n: int, m: int):
        """Braids the two anyons at positions 'n' and 'm'.
        If n < m, they are braided in a clockwise direction,
        if n > m, they are braided in a counterclockwise direction.

        Parameters
        ----------
        n : int
            The 1st anyon's position.
        m : int
            The 2nd anyon's position.

        Raises
        ------
        Exception
            Exceptions are raised if a braiding is attempted after a
            measurement, or if trying to braid two non-adjacent anyons.
        ValueError
            Is raised if the parameters are not strictly positive integers,
            or if incorrect positions are passed.

        Returns
        -------
        AnyonicCircuit
            A reference to the same circuit.

        """
        if self.__measured:
            raise Exception("System already measured! Cannot perform further braiding!")

        if not isinstance(m, int) or not isinstance(n, int):
            raise ValueError("n, m must be integers")

        if m < 1 or n < 1:
            raise ValueError("n, m must be higher than 0!")

        if m > self.__nb_anyons or n > self.__nb_anyons:
            ## Check if correct
            raise ValueError(
                f"The system has only {self.__nb_anyons} anyons! n, m are erroneous!"
            )

        if abs(n - m) != 1:
            raise Exception("You can only braid adjacent anyons!")

        if n < m:
            self.__unitary = self.__sigmas[n - 1] @ self.__unitary
        else:
            self.__unitary = self.__sigmas[m - 1].T.conjugate() @ self.__unitary

        self.__braids_history.append((n, m))

        self.__nb_braids += 1

        self.drawer.braid(m, n)

        return self

    def braid_sequence(self, braid: Sequence[Sequence[int]]):
        """Takes a sequence of [sigma operator, power], and applies the
        successive operators to the 'power'.
        The first operator in the sequence is the first to be applied.

        Parameters
        ----------
        braid : Sequence[Sequence[int]]
            A sequence of pairs of integers representing a braiding operator
            and an exponent.

        Raises
        ------
        ValueError
            Is raised if one of the operators' indices is not an integer
            greater or equal to 1, or is an incorrect index.

        Returns
        -------
        AnyonicCircuit
            A reference to the same circuit.
        """
        for ind, power in braid:
            if not isinstance(ind, int) or not isinstance(power, int):
                raise ValueError("Indices and powers must be integers!")
            if ind < 1:
                raise ValueError(
                    f"sigma_{ind} is not a valid braiding operator! "
                    f"The operators indices must be >= 1."
                )
            if ind >= self.__nb_anyons:
                raise ValueError(
                    f"sigma_{ind} is not a valid braiding operator! "
                    f"The operators indices must be < {self.__nb_anyons}."
                )
            # Computing m and n
            m = n = 0
            if power > 0:
                n = ind
                m = ind + 1
            elif power < 0:
                m = ind
                n = ind + 1
            else:  # if power=0, do nothing (identity)
                continue

            for _ in range(abs(power)):
                self.braid(n, m)

        return self

    def measure(self):
        """Performs a measurement on the whole circuit.

        Raises
        ------
        Exception
            Is raised if a measurement has already been carried.

        Returns
        -------
        AnyonicCircuit
            A reference to the same circuit.
        """
        if self.__measured:
            raise Exception("Cannot carry the measurements twice!")

        self.__measured = True
        self.drawer.measure()
        return self

    def history(self, output: str = "raw"):
        """Returns the history of all braiding operations that were performed
        in the circuit.
        Its output can either be the raw braiding operations (n, m), a list of
        braiding operators (sigmas), or a LaTeX string containing the product
        of all the braiding operators.

        Parameters
        ----------
        output : str, optional
            Can either be "raw", "sigmas", or "latex". The default is "raw".

        Raises
        ------
        ValueError
            Is raised if an incorrect output format is chosen.

        Returns
        -------
        List or String
            Either a list of (n, m) operations, a list of braiding operators,
            or a LaTeX string.

        """
        if not output in ["raw", "sigmas", "latex"]:
            raise ValueError('Output should be either: "raw", "sigmas" or "latex"')

        if output == "raw":
            return self.__braids_history

        elif output == "sigmas":
            ret = []
            for (n, m) in self.__braids_history:
                if m < n:
                    ret.append(f"is{m}")
                else:
                    ret.append(f"s{n}")
            return ret

        else:
            sigmas = self.history(output="sigmas")

            # Converting to a sigma notation with powers
            power_sigmas = []
            last_sigma = sigmas[-1] if len(sigmas) else None
            power = 0
            for sigma in reversed(sigmas):
                if sigma == last_sigma:
                    power += 1
                    continue
                # Done counting the powers of the last sigma, adding it
                power_sigmas.append((last_sigma, power))
                # resetting
                power = 1
                last_sigma = sigma

            # Handling the last sigma
            if power != 0:
                power_sigmas.append((last_sigma, power))

            # Converting to LaTeX
            latex = ""
            for sigma, p in power_sigmas:
                # Inverses of sigmas (negative powers)
                if sigma[0] == "i":
                    latex += r"\sigma_{" f"{sigma[2:]}" "}^{" f"{-p}" "}"
                # Sigmas (positive powers)
                else:
                    latex += r"\sigma_{" f"{sigma[1:]}" "}"
                    if p > 1:  # Only add exponents != 1
                        latex += "^{" f"{p}" "}"

            latex = "$ " + latex + "$"
            return latex

    def draw(self):
        """Draws the topological quantum circuit.

        Returns
        -------
        Figure
            The braid describing the topological quantum circuit.
        """
        return self.drawer.draw()

    def statevector(self) -> np.ndarray:
        """Computes and returns the current state vector of the circuit.

        Returns
        -------
        ndarray
            The state vector of the circuit.

        """
        return self.__unitary @ self.__initial_state

    def unitary(self) -> np.ndarray:
        """Returns the unitary representation of the quantum circuit.

        Returns
        -------
        ndarray
            The unitary matrix for the circuit.
        """
        return self.__unitary

    def run(self, shots: int = 1000):
        """Simulates the quantum circuit for 'shots' number of times and
        returns the measurement results.

        Parameters
        ----------
        shots : int, optional
            Number of times the circuit is simulated. The default is 1000.

        Raises
        ------
        Exception
            Is raised if the circuit is run without a measurement.

        Returns
        -------
        dict
            Contains the number of measurements for each measured state.

        """
        # Needs to be measured?
        if not self.__measured:
            raise Exception("The system was not measured!")

        statevector = self.statevector()
        probs = np.ravel(np.real(statevector * statevector.conjugate()))
        memory = np.random.choice(np.arange(self.__dim), p=probs, size=shots)

        idx, counts = np.unique(memory, return_counts=True)

        counts_dict = {}
        for i in range(len(idx)):
            counts_dict[str(idx[i])] = counts[i]

        return {"counts": counts_dict, "memory": memory}
