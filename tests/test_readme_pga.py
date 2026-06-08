from ase.build import molecule
from pointgroup import PointGroupAnalyzer, IdealGasThermoAuto

def test_pga_example():
    atoms = molecule("H2O")

    pga = PointGroupAnalyzer(atoms)

    assert pga.pointgroup is not None
    assert pga.symmetry_number > 0

def test_igta_example():
    atoms = molecule("H2")

    vibrations = [0.5]

    igta = IdealGasThermoAuto(atoms, vibrations,
                              potentialenergy=0)

    G = igta.get_gibbs_energy(
        temperature=298.15,
        pressure=101325,
    )

    assert isinstance(G, float)