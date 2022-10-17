from typing import List, Tuple, Sequence
import numpy as np

from .lib.basis_generator import gen_basis
from .lib.get_generators import braiding_generator
from .lib.drawer import Drawer


class AnyonicCircuit:
    def __init__(self, nb_qudits: int, nb_anyons_per_qudit: int):
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
        return self.__nb_qudits

    @property
    def nb_anyons_per_qudits(self):
        return self.__nb_anyons_per_qudit

    @property
    def drawer(self):
        return self.__drawer

    @property
    def dim(self):
        return self.__dim

    @property
    def basis(self):
        return self.__basis

    def __get_basis(self) -> Tuple[np.ndarray, int]:
        basis = gen_basis(self.__nb_qudits, self.__nb_anyons_per_qudit)

        return basis, len(basis)

    def __get_sigmas(self) -> List[np.ndarray]:
        sigmas = []
        for index in range(1, self.__nb_anyons):
            sigma = braiding_generator(
                index, self.__nb_qudits, self.__nb_anyons_per_qudit
            )
            sigmas.append(np.array(sigma))
        return sigmas

    def initialize(self, input_state: np.ndarray):
        if self.__nb_braids > 0:
            raise Exception(
                "Initialization should happen before any braiding operation is performed!"
            )

        input_state = np.array(input_state)
        if not np.size(input_state) == self.__dim:
            raise ValueError("The state has wrong dimension. Should be {self.__dim}")

        norm = np.sum(np.real(input_state * input_state.conjugate()))
        if not np.isclose(norm, 1, 5):
            raise ValueError("The input state is not normalized correctly!")

        self.__initial_state = np.reshape(input_state, (self.__dim, 1))

        return self

    def braid(self, n: int, m: int):
        if self.__measured:
            raise Exception("System already measured! Cannot perform further braiding!")

        if abs(n - m) != 1:
            raise Exception("You can only braid adjacent anyons!")

        if not isinstance(m, int) or not isinstance(n, int):
            raise ValueError("n, m must be integers")
        if m < 1 or n < 1:
            raise ValueError("n, m must be higher than 0!")

        if m > self.__nb_anyons or n > self.__nb_anyons:
            ## Check if correct
            raise ValueError(
                f"The system has only {self.__nb_anyons} anyons! n, m are erroneous!"
            )

        if n < m:
            self.__unitary = self.__sigmas[n - 1] @ self.__unitary
        else:
            self.__unitary = self.__sigmas[m - 1].T.conjugate() @ self.__unitary


        self.__braids_history.append((m, n))

        self.__nb_braids += 1

        self.drawer.braid(m, n)

        return self

    def braid_sequence(self, braid: Sequence[Sequence[int]]):
        """ Takes a sequence of [sigma operator, power], and applies the
        successive operators to the 'power'.
        The first operator in the sequence is the first to be applied.
        """
        for ind, power in braid:
            if not isinstance(ind, int) or not isinstance(power, int):
                raise ValueError("Indices and powers must be integers!")
            if ind < 1:
                raise ValueError(f"sigma_{ind} is not a valid braiding operator! "
                                 f"The operators indices must be >= 1.")
            if ind >= self.__nb_anyons:
                raise ValueError(f"sigma_{ind} is not a valid braiding operator! "
                                 f"The operators indices must be < {self.__nb_anyons}.")
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
        self.__measured = True
        self.drawer.measure()
        return self

    def history(self, output="raw"):
        if not output in ["raw", "sigmas", "latex"]:
            raise ValueError('Output should be either: "raw", "sigmas" or "latex"')

        if output == "raw":
            return self.__braids_history

        elif output == "sigmas":
            ret = []
            for (m, n) in self.__braids_history:
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
                    latex += "\sigma_{"f"{sigma[2:]}""}^{"f"{-p}""}"
                # Sigmas (positive powers)
                else:
                    latex += "\sigma_{"f"{sigma[1:]}""}"
                    if p > 1:  # Only add exponents != 1
                        latex += "^{"f"{p}""}"

            latex = "$ " + latex + "$"
            return latex

    def draw(self):
        return self.drawer.draw()

    def statevector(self) -> np.ndarray:
        return self.__unitary @ self.__initial_state

    def unitary(self) -> np.ndarray:
        return self.__unitary

    def run(self, shots=1024):
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
