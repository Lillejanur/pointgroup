from ase.thermochemistry import IdealGasThermo, AbstractMode
from ase import Atoms
from ase.units import invcm
from collections.abc import Sequence
from typing import Literal
import warnings
import numpy as np
from .pointgroup import PointGroupAnalyzer

_VIB_SELECT_OPTIONS = Literal['exact', 'all', 'highest', 'abs_highest']
_TRIPLET_GS_MOLECULES = {'O2', 'CH2'} # If expanded, check the detection code

class IdealGasThermoAuto(IdealGasThermo):
    """
    Create an IdealGasThermo object with automatic symmetry analysis.

    The molecular geometry and symmetry number are determined from
    the atomic structure using PointGroupAnalyzer.

    Inputs:

    atoms : an ASE Atoms object
        used to calculate symmetries, rotational moments of inertia and
        molecular mass, as well as to retrieve energy if potentialenergy is
        unspecified and to check spin if unspecified.
    vibrations : sequence of complex
        a list of the vibrational frequencies of the molecule; default unit is
        eV. See vib_energies in ASE's IdealGasThermo for details.
    vib_unit : 'eV' (default), 'invcm', 'cm^-1'
        unit for vibrational frequencies. 'invcm' and 'cm^-1' are identical.
    fmax : float or None
        maximum force; defaults to 0.1 eV/Å. If atoms have forces
        available and fmax is exceeded, an error will be raised. The point
        is to catch unfinished relaxations.
    potentialenergy : float or None
        the potential energy in eV. If None (default), a cached energy from the
        Atoms' calculator will be used.
    vib_selection : 'exact', 'highest' (default), 'abs_highest', 'all'
        selection of input vibrational energies considered to be true
        vibrations (excluding translations and rotations) implied by the
        geometry and number of atoms. See IdealGasThermo for details.
    ignore_imag_modes : bool
        see IdealGasThermo for details.
    modes : AbstractMode
        see IdealGasThermo for details (seems unused).
    spin : float or None
        the total electronic spin; 0 for molecules in which all
        electrons are paired, 0.5 for a free radical with a single unpaired
        electron, 1.0 for a triplet with two unpaired electrons. If None,
        it will be interpreted as 0 and atoms will be checked if non-zero spin
        is reasonable and if so raise an error (which can be overridden with
        spin=0).
    eigtol : float
        tolerance of inertia eigenvalues, normalized by trace. Half the
        tolerance of other codes like pymatgen or pypi pointgroup, as
        those normalize by half the trace.
    angtol : float
        angle tolerance in degrees.
    disttol : float
        distance tolerance.
    hardtol : float
        other tolerances.
    """

    def __init__(self, atoms: Atoms,
                 vibrations: Sequence[complex],
                 vib_unit: Literal['eV', 'invcm', 'cm^-1'] = 'eV',
                 fmax: float = 0.1,
                 potentialenergy: float | None = None,
                 spin: float | None = None,
                 vib_selection: _VIB_SELECT_OPTIONS | None = 'highest',
                 ignore_imag_modes: bool = False,
                 modes: Sequence[AbstractMode] | None = None,
                 eigtol: float = 0.005,
                 angtol: float = 4.0,
                 disttol: float = 0.2,
                 hardtol: float = 1e-6) -> None:

        # We have to use the input atoms, rather than a copy, because
        # the copy does not transfer the calculator
        calc = getattr(atoms, 'calc', None)
        results = getattr(calc, 'results', {}) if calc is not None else {}
        # Use cached energy if it exists, otherwise 0
        if potentialenergy is None:        
            potentialenergy = results.get('energy')
            if potentialenergy is None:
                raise ValueError('No cached energy found in atoms. '
                                 'Please specify potentialenergy explicitly.'
                                )
        # Check that existing forces do not exceed fmax
        if fmax is not None:
            forces = results.get("forces")
            if forces is not None:
                max_force = np.linalg.norm(forces, axis=1).max()
                if max_force > fmax:
                    raise ValueError(
                        f'fmax={fmax} eV/Å exceeded. '
                        f'(max force = {max_force:.3f} eV/Å). '
                        'Atoms appear not to be relaxed.'
                    )
        self.atoms = atoms.copy()
        if spin is None:
            num_electrons = self.atoms.get_atomic_numbers().sum()
            # Formula may be unstable for expansion of triplet ground state
            # molecules.
            formula = atoms.get_chemical_formula()
            if num_electrons % 2:
                raise ValueError(
                    'Odd number of electrons; at least one unpaired. '
                    'Please specify spin explicitly.'
                )
            if formula in _TRIPLET_GS_MOLECULES:
                raise ValueError(
                    f'{formula} is commonly in a triplet ground-state. '
                    'Please specify spin explicitly.'
                )
            spin = 0
        elif spin < 0:
            raise ValueError('Total electronic spin must be non-negative.')

        # Appearently, the ASE class do not accept periodic boundary conditions
        self.atoms.set_pbc(False)
        pga = PointGroupAnalyzer(self.atoms, eigtol=eigtol,
                                 angtol=angtol, disttol=disttol,
                                 hardtol=hardtol)
        self.pga = pga
        
        if vib_unit in ['invcm', 'cm^-1']:
            if np.max(np.abs(vibrations)) < 1:
                warnings.warn(
                    'Largest frequency below 10 cm^-1 is unusual. '
                    'Did you provide frequencies in eV?'
                )
            vib_energies = [v * invcm for v in vibrations]
        elif vib_unit == 'eV':
            if np.max(np.abs(vibrations)) > 10:
                warnings.warn(
                    'Vibrational energies above 10 eV are unusual. '
                    'Did you provide frequencies in cm^-1?',
                    UserWarning
                )
            vib_energies = vibrations
        else:
            raise ValueError(
                f'Unknown vib_unit={vib_unit!r}. '
                "Expected 'eV', 'invcm', or 'cm^-1'.",
                UserWarning
            )

        super().__init__(
            vib_energies=vib_energies,
            geometry=pga.geometry,
            potentialenergy=potentialenergy,
            atoms=self.atoms,
            symmetrynumber=pga.symmetry_number,
            spin=spin,
            vib_selection=vib_selection,
            ignore_imag_modes=ignore_imag_modes,
            modes=modes
        )

    #def __getattr__(self, name):
    #    return getattr(self.idealgasthermo, name)