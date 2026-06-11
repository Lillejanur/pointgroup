# PointGroup

A Python package for automatic molecular point group detection and ideal-gas thermochemistry built on top of the Atomic Simulation Environment (ASE).

It combines symmetry analysis with ASE’s thermochemistry tools to reduce manual setup and prevents errors in molecular thermochemistry workflows.

---

## Features

- Automatic point group detection from atomic geometry
- Symmetry number assignment
- Geometry classification (monoatomic / linear / nonlinear)
- Wrapper around ASE `IdealGasThermo`

---

## Requirements

- Python ≥ 3.10
- numpy
- scipy
- ase

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/Lillejanur/pointgroup.git
cd pointgroup
pip install -e .
```
---

## Example usage point group detection

```
python

from ase.build import molecule
from pointgroup import PointGroupAnalyzer

atoms = molecule("H2O")

pga = PointGroupAnalyzer(atoms)

print(pga.pointgroup) # C2v
print(pga.symmetry_number) # 2
print(pga.geometry) # 'linear'
```

## Example usage ideal gas thermo wrapper

```
python

from ase.build import molecule
from pointgroup import IdealGasThermoAuto

atoms = molecule("H2")

vib_energies = [0.5]  # eV; cm^-1 works with vib_unit='invcm' or vib_unit='cm^-1'

thermo = IdealGasThermoAuto(atoms, vib_energies, potentialenergy=0)

G = thermo.get_gibbs_energy(
    temperature=298.15,
    pressure=101325
)
```

