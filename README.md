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

### Basic Example
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


## License

Copyright © 2022, [Constantine Quantum Technologies](https://cqtech.org). Released under the [Apache License 2.0](LICENSE).
