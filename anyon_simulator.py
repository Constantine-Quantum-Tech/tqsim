import numpy as np
from IPython.display import display, Latex
import matplotlib.pyplot as plt
import matplotlib as mpl
import braiding_generators.fib_multi_qudits_ as fmq
from copy import deepcopy

class AnyonSimulator():

    def __init__(self, n_qubits=1, n_anyons_per_qubit=3):

        if n_anyons_per_qubit == 3:
            if n_qubits > 4:
                raise ValueError("Number of qubits is limited to 4")
        elif n_anyons_per_qubit == 4:
            if n_qubits > 3:
                raise ValueError("Number of qubits is limited to 3")
        else:
            raise ValueError("Number of anyons per qubit can be either 3 or 4")

        # Number of qubits and anyons
        self.n_qubits = n_qubits
        self.n_anyons_per_qubit = n_anyons_per_qubit
        self.nb_anyons = n_anyons_per_qubit * n_qubits

        # Current braiding ops in the circuit
        self.braids = 0
        self._braids_history = []

        # Measurement flag
        self._measured = False

        # Braids drawing data
        self._i = np.linspace(0, 1)
        self._il = np.linspace(0, 0.46, 25)
        self._ih = np.linspace(0.54, 1, 25)
        self._id = self._ih[0] - self._il[-1]
        self._x = dict()
        self._y = dict()
        self.labels = {}

        for i in range(self.n_qubits):
            for j in range(self.n_anyons_per_qubit):
                top_current_qubit = (n_anyons_per_qubit + 1) * n_qubits - i * (n_anyons_per_qubit + 1)
                self._x[f"{i + 1}_{j+1}"] = [np.copy(self._i)]
                self._y[f"{i + 1}_{j+1}"] = [ (top_current_qubit - j) * np.ones_like(self._i)]
                self.labels[f"{i + 1}_{j+1}"] = {'x': self._x[f"{i + 1}_{j+1}"][0][0]-0.75,
                                                 'y': self._y[f"{i + 1}_{j+1}"][0][0]}
        # Unused?
        # self._omega = np.exp(1j * np.pi / 5)
        # self._phi = (1 + np.sqrt(5)) / 2

        # Braiding operators
        self._s = []
        for ii in range(1, self.nb_anyons):
            self._s.append(np.array(fmq.braiding_generator(index=ii,
                                                            n_qudits=self.n_qubits,
                                                            qudit_len=self.n_anyons_per_qubit-1,
                                                            show=False)))
        # Basis of the Hilbert space and its dimension
        self._basis = fmq.find_basis(n_qudits=self.n_qubits,
                                      qudit_len=self.n_anyons_per_qubit-1)
        self._dim = len(self._basis)

        # Initial unitary of the circuit (identity)
        self._unitary = np.eye(self._dim)

        self._anyon_pos = [p for p in range(self.nb_anyons)]

        #print(f'Fusion space is {self._dim} dimensional.')

    def _sigmoid(self, x):
        return 1 / (1 + np.exp(-(x - 0.5) * 12))

    def braid(self, n, m):
        if self._measured:
            raise Exception("Cannot braid after making a measurement!")

        k = [k for k in self._x.keys()]

        qubit1 = (n - 1) // self.n_anyons_per_qubit + 1
        qubit2 = (m - 1) // self.n_anyons_per_qubit + 1
        qn = (n - 1) % self.n_anyons_per_qubit + 1
        qm = (m - 1) % self.n_anyons_per_qubit + 1

        top = f"{qubit1}_{qn}"
        bottom = f"{qubit2}_{qm}"



        d = abs(self._y[top][-1][-1] - self._y[bottom][-1][-1])

        if abs(n - m) != 1 or not (top in k and bottom in k):
            raise ValueError("Wrong n, m!")

        k.remove(top)
        k.remove(bottom)

        # k idle(s)
        for i in k:
            self._x[i].append(self._x[i][-1][-1] + self._i)
            self._y[i].append(self._y[i][-1][-1] * np.ones_like(self._i))


#         # n above
        self._x[top].append(self._x[top][-1][-1] + self._i)
        self._y[top].append(self._y[top][-1][-1] + d * np.sign(n - m) * self._sigmoid(self._i))

#         # m below

        b_y = self._y[bottom][-1][-1]
        a_xl = self._x[bottom][-1][-1] + self._il
        a_yl = self._y[bottom][-1][-1] - d * np.sign(n - m) * self._sigmoid(self._il)

        self._x[bottom].append(a_xl)
        self._y[bottom].append(a_yl)

        a_xh = self._x[bottom][-1][-1] + self._id + self._il
        a_yh = b_y - d * np.sign(n - m) * self._sigmoid(self._ih)

        self._x[bottom].append(a_xh)
        self._y[bottom].append(a_yh)

        self._x[top], self._x[bottom] = self._x[bottom], self._x[top]
        self._y[top], self._y[bottom] = self._y[bottom], self._y[top]

        self._anyon_pos[n-1], self._anyon_pos[m-1] = self._anyon_pos[m-1], self._anyon_pos[n-1]

        self.braids += 1

        if m == n - 1:
            self._unitary = self._s[m-1].T.conjugate() @ self._unitary ####
            self._braids_history.append(f"is{m}")
        else:
            self._unitary = self._s[n-1] @ self._unitary ####
            self._braids_history.append(f"s{n}")


    def _identity(self):
        for i in self._x.keys():
            self._x[i].append(self._x[i][-1][-1] * np.ones_like(self._i))
            self._y[i].append(self._y[i][-1][-1] * np.ones_like(self._i))


    def measure(self):
        self._identity()

        for i in range(self.n_qubits):
            for j in range(self.n_anyons_per_qubit - 1):
                # Identities
                for k in range(j + 3, self.n_anyons_per_qubit + 1):
                    self._x[f"{i + 1}_{k}"].append(self._x[f"{i + 1}_{k}"][-1][-1] + self._i)
                    self._y[f"{i + 1}_{k}"].append(self._y[f"{i + 1}_{k}"][-1][-1] * np.ones_like(self._i))
                # Fusion
                self._x[f"{i + 1}_{j+1}"].append(self._x[f"{i + 1}_{j+1}"][-1][-1] + self._i)
                self._x[f"{i + 1}_{j+2}"].append(self._x[f"{i + 1}_{j+2}"][-1][-1] + self._i)
                d = 0.5 * (self._y[f"{i + 1}_{j+1}"][-1][-1] - self._y[f"{i + 1}_{j+2}"][-1][-1])
                self._y[f"{i + 1}_{j+1}"].append(self._y[f"{i + 1}_{j+1}"][-1][-1]  - d * self._sigmoid(self._i))
                self._y[f"{i + 1}_{j+2}"].append(self._y[f"{i + 1}_{j+2}"][-1][-1]  + d * self._sigmoid(self._i))

        self._measured = True

    def draw_circuit(self, colors=dict()):
        top_current_qubit = (self.n_anyons_per_qubit + 1) * self.n_qubits

        new_colors = {}

        for k in colors.keys():
            new_colors[k] = k[:-1] + str(self._anyon_pos.index(int(k[-1])-1) + 1)


        ccolors = {}
        for key in colors.keys():
            ccolors[new_colors[key]] = colors[key]

        for k in self._x.keys():
            if ccolors.get(k) is None:
                ccolors[k] = "black"


        # if self.braids > 0:
            # width, height = plt.figaspect(1 / self.braids)
            # fig = plt.figure(figsize=(width, height))

        mpl.rcParams['figure.figsize'] = (1 * self.braids, 6)

        plt.ylim([0, top_current_qubit + 1])

        for b in self._x.keys():
            # Anyons' worldlines
            # Anyons' labels
            plt.text(self.labels[b]['x'], self.labels[b]['y'], b)
            for i, _ in enumerate(self._x[b]):
                plt.plot(self._x[b][i], self._y[b][i], ccolors[b], linewidth=2)
        plt.axis("off")
        plt.show()

    def get_braids_history(self):
        latex = " ".join(self._braids_history[::-1])
        latex = latex.replace("is", "\sigma^{-1}_")
        for ii in range(1, self.nb_anyons):
            latex = latex.replace(f"s{ii}", f"\sigma_{ii}")
        display(Latex(f'${latex}$'))

    def statevector(self, input_state):

        norm = 0
        for amplitude in input_state:
            norm += amplitude * amplitude.conjugate()

        if not np.isclose(norm, 1, 5):
            raise 'Wrong state!'

        output_statevector = self._unitary @ input_state

        return output_statevector

    def unitary(self):
        fmq.cplot(self._unitary)
        return self._unitary

    def simulate(self, input_state, shots=1):
        amplitudes = self.statevector(input_state)

        probs = []
        for amp in amplitudes:
            probs.append((amp * amp.conjugate()).real)
        return np.random.choice([ii for ii in range(self._dim)], p=probs, size=shots)
