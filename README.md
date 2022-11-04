<img align="right" width="70" src="https://raw.githubusercontent.com/Constantine-Quantum-Tech/tqsim/main/images/cqtech_logo.png" alt="CQTech">

# TQSim

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/Constantine-Quantum-Tech/tqsim/tree/main.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/Constantine-Quantum-Tech/tqsim/tree/main)
[![Version](https://img.shields.io/pypi/v/tqsim?style=flat-square)](https://pypi.org/project/tqsim/)
[![License](https://img.shields.io/pypi/dm/tqsim?style=flat-square)](https://pypi.org/project/tqsim/)
[![License](https://img.shields.io/github/license/Constantine-Quantum-Tech/tqsim?style=flat-square)](LICENSE)


TQSim stands for Topological Quantum Simulator. It is an open-source library for simulating topological quantum computer based on anyons. 


## Installation

You can install TQSim from pip using

```bash
pip install --upgrade tqsim
```

## Usage

### 1. Basic Example
In this example, we create a circuit with 2 qudits, made of 3 anyons each. We then braid the anyons manually.
```python
from tqsim import AnyonicCircuit

circuit = AnyonicCircuit(2, 3) # Create a circuit having 2 qudits and 3 anyons per qudits
circuit.braid(1, 2) # Braids the first with the second anyon
circuit.braid(3, 4) # Braids the first with the second anyon
circuit.braid(2, 1)
circuit.measure() # Measure the system by fusing the anyons
circuit.draw() # Draw the circuit
```
Here is the output of the `draw()` method:

![Circuit Output](https://i.ibb.co/3z9pFmQ/example.png)

Simulating the circuit:

```python
circuit.run(shots = 50)

```

Output:
```bash
{'counts': {'0': 16, '2': 20, '4': 14}, 'memory': array([4, 2, 2, 2, 2, 2, 2, 4, 2, 4, 2, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4,
       2, 4, 0, 2, 0, 0, 0, 0, 4, 4, 2, 2, 2, 4, 2, 2, 0, 0, 2, 4, 2, 2,
       4, 2, 4, 4, 0, 2])}
```

### 2. Simulating a Hadamard gate
Here we simulate the application of a Hadamard gate on a single qudit with 3 anyons.
Unlike the previous example, we will use a braiding sequence of braiding operators and their corresponding powers.
```python
from tqsim import AnyonicCircuit

circuit = AnyonicCircuit(nb_qudits=1, nb_anyons_per_qudit=3)  # Create a circuit with 1 qudit composed of 3 anyons
circuit.initialize([0,0,1])  # We initialize the circuit in the last state (state 2).
                            # For this circuit, we have 3 basis states: [0, 1, 2].

# The Hadamard gate braiding sequence in terms of braiding operators
had_sequence = [[1, 2], [2, 2], [1, -2], [2, -2], [1, 2], [2, 4], [1, -2], [2, 2],
                [1, 2], [2, -2], [1, 2], [2, -2], [1, 4]]

circuit.braid_sequence(had_sequence)  # We apply the braiding sequence.
                                      # This should put our qudit in a superposition of the states 2 and 1.
circuit.measure()  # Measure the system by fusing the anyons
circuit.draw()  # Draw the circuit
```
The Hadamard braid looks like this

![Circuit Output](https://i.ibb.co/5kDVWgf/example-hadamard.png)

Simulating the circuit:

```python
result = circuit.run(shots = 1000)  # Run the circuit 1000 times.
print(result['counts'])  # Show the counts only.

```

Output:
```bash
{'1': 493, '2': 507}
```

## Authors and Citation
*Abdellah Tounsi, Mohamed Messaoud Louamri, Nacer eddine Belaloui, and Mohamed Taha Rouabah – Constantine Quantum Technologies.*

If you have used TQSim in your work, please use the [BibTeX file](citation.bib) to cite the authors.

## License

Copyright © 2022, [Constantine Quantum Technologies](https://cqtech.org). Released under the [Apache License 2.0](LICENSE).
