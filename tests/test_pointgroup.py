from dataclasses import dataclass

import pytest

from pathlib import Path
from ase.build import molecule
from ase.io import read
from pointgroup import PointGroupAnalyzer

@dataclass
class MoleculeData:
    label: str
    pointgroup: str
    symmetry_number: int

    @property
    def filename(self):
        return f'{self.label}.xyz'

ase_database_cases = [
    MoleculeData('Al', 'Kh', 1),
    MoleculeData('CH3CONH2', 'C1', 1),
    MoleculeData('HCOOH', 'Cs', 1),
    MoleculeData('HOCl', 'Cs', 1),
    MoleculeData('CH3CHO', 'Cs', 1),
    MoleculeData('CO', 'C*v', 1),
    MoleculeData('N2H4', 'C2', 2),
    MoleculeData('H2O', 'C2v', 2),
    # C4H4NH and COF2 can be mistakenly identified as symmetric
    # due to mass distribution, which might lead to incorrect
    # identification as point group Cs
    MoleculeData('C4H4NH', 'C2v', 2),
    MoleculeData('COF2', 'C2v', 2),
    MoleculeData('butadiene', 'C2h', 2),
    MoleculeData('NH3', 'C3v', 3),
    MoleculeData('C2H4', 'D2h', 4),
    MoleculeData('C5H8', 'D2d', 4),
    MoleculeData('C2H6', 'D3d', 6),
    MoleculeData('C6H6', 'D6h', 12),
    MoleculeData('CH4', 'Td', 12),
    MoleculeData('C60', 'Ih', 60),
]

@pytest.mark.parametrize(
    'moldata',
    ase_database_cases,
    ids=[case.label for case in ase_database_cases]
)
def test_pointgroup(moldata, rotate=[30, 30, 30]):
    label = moldata.label
    try:
        mol = molecule(label)
    except Exception as e:
        pytest.fail(f"ASE does not recognize molecule label '{label}': {e}")
    # Rotate to break any accidental symmetries
    mol.euler_rotate(*rotate)
    pga = PointGroupAnalyzer(mol)
    pg_calc = pga.pointgroup
    sym_calc = pga.symmetry_number
    assert pg_calc == moldata.pointgroup
    assert sym_calc == moldata.symmetry_number

# --------------Molecules from structure files in testdata---------------
TESTDATA = Path(__file__).parent / "testdata"

def load_molecule(moldata, rotate=[30, 30, 30]):
    mol = read(TESTDATA / moldata.filename)
    mol.set_pbc(False)

    if rotate:
        mol.euler_rotate(*rotate)

    return mol

cases_from_testdata = [
    MoleculeData(label='SF6', pointgroup='Oh', symmetry_number=24),
    MoleculeData(label='B12H12', pointgroup='Ih', symmetry_number=60),
    # HOBr can be incorrectly detected as linear because of its mass
    # distribution
    MoleculeData(label='HOBr', pointgroup='Cs', symmetry_number=1),
    MoleculeData(label='cyanoacetylene', pointgroup='C*v', symmetry_number=1),
    MoleculeData(label='C4H2', pointgroup='D*h', symmetry_number=2),
    MoleculeData(label='C3O2', pointgroup='D*h', symmetry_number=2),
    MoleculeData(label='boric_acid', pointgroup='C3h', symmetry_number=3),
    MoleculeData(label='twistane', pointgroup='D2', symmetry_number=4),
    MoleculeData(label='XeF4O', pointgroup='C4v', symmetry_number=4),
    MoleculeData(label='C60F36', pointgroup='T', symmetry_number=12),
    MoleculeData(label='1,2-dichloro-1,2-difluoroethane',
                 pointgroup='Ci', symmetry_number=1),
    MoleculeData(label='cubane', pointgroup='Oh', symmetry_number=24),
    MoleculeData(label='C5H4F4', pointgroup='S4', symmetry_number=2),
    MoleculeData(label='uranocene', pointgroup='D8h', symmetry_number=16),
    MoleculeData(label='S8', pointgroup='D4d', symmetry_number=8),
    MoleculeData(label='XeF4', pointgroup='D4h', symmetry_number=8),
    MoleculeData(label='corannulene', pointgroup='C5v', symmetry_number=5),
    MoleculeData(label='ferrocene', pointgroup='D5d', symmetry_number=10),
    MoleculeData(label='C70', pointgroup='D5h', symmetry_number=10),
]


@pytest.mark.parametrize(
    'moldata',
    cases_from_testdata,
    ids=[case.label for case in cases_from_testdata],
)
def test_pointgroups_in_testdata(moldata):
    mol = load_molecule(moldata)
    pga = PointGroupAnalyzer(mol)
    assert pga.pointgroup == moldata.pointgroup
    assert pga.symmetry_number == moldata.symmetry_number


# Some molecules with slight deviations from perfect symmetry
imperfect_molecules = [
    MoleculeData(label='C20', pointgroup='Ih', symmetry_number=60),
    MoleculeData(
        label='thorium_nitrate_ion', pointgroup='Th', symmetry_number=12
    ),
]


@pytest.mark.parametrize('moldata', imperfect_molecules)
def test_pointgroups_imperfect(moldata):
    mol = load_molecule(moldata)
    mol.euler_rotate(30.0, 30.0, 30.0)
    pga = PointGroupAnalyzer(mol, eigtol=0.015, angtol=6, disttol=0.3)
    assert pga.pointgroup == moldata.pointgroup
    assert pga.symmetry_number == moldata.symmetry_number
