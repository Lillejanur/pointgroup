import pytest
from pathlib import Path
from ase.build import molecule
from ase.io import read
from ase.units import invcm
from ase.thermochemistry import IdealGasThermo
from pointgroup import IdealGasThermoAuto

T = 300
p = 1e5

moldata = [
    ('H2', [4400], 'linear', 2),
    ('CO', [2100], 'linear', 1),
    ('H2O', [1600, 3700, 3800], 'nonlinear', 2),
    ('CH4', [1280, 1290, 1300, 1510, 1520, 2900, 3000, 3020, 3025],
     'nonlinear', 12),
]

@pytest.mark.parametrize(
    'label, freqs, geom, sym',
    moldata,
    ids=[mol[0] for mol in moldata]
)
def test_mols_from_table(label, freqs, geom, sym, capsys):
    mol = molecule(label)
    vibs = [freq * invcm for freq in freqs]
    igta = IdealGasThermoAuto(mol, vibs, potentialenergy=0)

    assert igta.pga.geometry == geom
    assert igta.pga.symmetry_number == sym

    G_auto = igta.get_gibbs_energy(temperature=T, pressure=p)
    out_auto = capsys.readouterr().out

    igt = IdealGasThermo(vibs, geom, atoms=mol,
                         symmetrynumber=sym, spin=0)
    G_ref = igt.get_gibbs_energy(temperature=T, pressure=p)
    out_ref = capsys.readouterr().out
    assert G_auto == pytest.approx(G_ref, rel=1e-12, abs=1e-12)
    assert out_auto == out_ref # May be fragile

@pytest.mark.parametrize(
    'label, freqs, geom, sym',
    moldata,
    ids=[mol[0] for mol in moldata]
)
def test_freq_handling(label, freqs, geom, sym):
    # Enter frequencies in invcm, program expects eV
    mol = molecule(label)
    high_vib_warning = r"Vibrational energies above 10 eV are unusual"
    with pytest.warns(UserWarning, match=high_vib_warning):
        IdealGasThermoAuto(mol, freqs, potentialenergy=0)

    # Enter frequencies in eV, program expects invcm
    vib_energies = [freq * invcm for freq in freqs]
    low_vib_warning = r"Largest frequency below 10 cm"
    with pytest.warns(UserWarning, match=low_vib_warning):
        IdealGasThermoAuto(mol, vib_energies, vib_unit='invcm',
                           potentialenergy=0)

    # Enter frequencies in invcm correctly. Compare with reference.

    igta3 = IdealGasThermoAuto(mol, freqs, vib_unit='cm^-1', potentialenergy=0)
    G3 = igta3.get_gibbs_energy(T, p)
    igt = IdealGasThermo(vib_energies,
                         geometry=geom,
                         atoms=mol,
                         symmetrynumber=sym,
                         spin=0)
    G_ref = igt.get_gibbs_energy(T, p)
    assert G3 == pytest.approx(G_ref, rel=1e-12, abs=1e-12)

moldata2 = [
    #('H2', [0.5], 0),
    ('O2', [0.3], 1, 'commonly in a triplet', 'linear', 2),
    ('CH2_s3B1d', [0.1, 0.2, 0.2], 1, 'commonly in a triplet', 'nonlinear', 2),
    ('CN', [0.35], 0.5, 'Odd number of electrons', 'linear', 1),
]

@pytest.mark.parametrize(
    'label, vib_energies, spin, match, geom, sym',
    moldata2,
    ids=[mol[0] for mol in moldata2]
)
def test_spin(label, vib_energies, spin, match, geom, sym):
    mol = molecule(label)
    with pytest.raises(ValueError, match=match):
        IdealGasThermoAuto(mol, vib_energies, potentialenergy=0)

    # Override with spin=0 should not raise errors
    IdealGasThermoAuto(mol, vib_energies, spin=0, potentialenergy=0)

    # Check that correct spin produces same result as IdealGasThermo
    igta = IdealGasThermoAuto(mol, vib_energies, spin=spin, potentialenergy=0)
    G = igta.get_gibbs_energy(T, p)
    igt = IdealGasThermo(vib_energies,
                         geometry=geom,
                         atoms=mol,
                         symmetrynumber=sym,
                         spin=spin)
    G_ref = igt.get_gibbs_energy(T, p)
    assert G == pytest.approx(G_ref, rel=1e-12, abs=1e-12)

TESTDATA = Path(__file__).parent / "testdata"
rlx_mols = ('H2', 'H2O', 'NH3')


@pytest.mark.parametrize(
    'label',
    rlx_mols,
    ids=rlx_mols
)
def test_cached_energy_forces(label):
    filename = label + '_rlx.traj'
    mol = read(TESTDATA / filename)
    freqs = mol.info['freqs']
    igta = IdealGasThermoAuto(mol, vibrations=freqs, vib_unit='invcm')
    G = igta.get_gibbs_energy(T, p)

    geom = igta.pga.geometry
    sym = igta.pga.symmetry_number
    energy = mol.calc.results['energy']

    # IdealGasThermo cannot accept periodic boundary conditions
    mol_nopbc = mol.copy()
    mol_nopbc.set_pbc(False)

    vib_energies = [freq * invcm for freq in freqs]
    igt = IdealGasThermo(vib_energies,
                         geometry=geom,
                         potentialenergy=energy,
                         atoms=mol_nopbc,
                         symmetrynumber=sym,
                         spin=0)
    G_ref = igt.get_gibbs_energy(T, p)
    assert G == pytest.approx(G_ref, rel=1e-12, abs=1e-12)

def test_too_high_force():
    # This structure has too high max force
    mol = read(TESTDATA / 'CH3OH_rlx.traj')
    freqs = mol.info['freqs']
    with pytest.raises(ValueError, match=r"max force ="):
        IdealGasThermoAuto(mol, vibrations=freqs, vib_unit='invcm')