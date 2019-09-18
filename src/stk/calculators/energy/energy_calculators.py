"""
Energy Calculators
==================

#. :class:`.MMFFEnergy`
#. :class:`.UFFEnergy`
#. :class:`.MacroModelEnergy`
#. :class:`.MOPACEnergy`
#. :class:`.FormationEnergy`
#. :class:`.XTBEnergy`
#. :class:`.XTBFreeEnergy`


Energy calculators are objects which calculate the energy of molecules.
Each :class:`.EnergyCalculator` is initialized with some settings and
calculates the energy of a molecule with
:meth:`~.EnergyCalculator.get_energy`.

.. code-block:: python

    import stk

    # Energy calculators work with any Molecule objects, such as
    # BuildingBlock or ConstructedMolecule.
    mol1 = stk.BuildingBlock('[Br]CC[Br]', ['bromine'])

    chain = stk.polymer.Linear('A', [0], 12)
    mol2 = stk.ConstructedMolecule([mol1], chain)

    # Create the energy calculator.
    mmff = stk.MMFFEnergy()

    # Calculate energies of various molecules.
    mol1_energy = mmff.get_energy(mol1)
    mol2_energy = mmff.get_energy(mol2)

By default, calling :meth:`~.EnergyCalculator.get_energy` twice on the
same molecule will calculate the energy a second time. However, we can
use the `use_cache` option to prevent
recalculations when the same molecule is given to the same energy
calculator a second time

.. code-block:: python

    caching_mmff = MMFFEnergy(use_cache=True)
    # Calculate the energy the first time.
    energy1 = caching_mmff.get_energy(mol1)
    # The second time, the energy is returned directly from memory, a
    # second calculation is not run.
    energy2 = caching_mmff.get_energy(mol1)


.. _`adding energy calculators`:

Making New Energy Calculators
-----------------------------

New energy calculators can be made by simply making a class which
inherits the :class:`.EnergyCalculator` class. In addition to this,
the new class must define a :meth:`~.EnergyCalculator.get_energy`
method. The method must take 1 parameter, `mol`. The method will then
calculate and return the energy. There are no requirements regarding
how it should go about calculating the energy. New energy calculators
can be added into the :mod:`.energy_calculators` submodule or into a
new submodule.

"""

import rdkit.Chem.AllChem as rdkit
import logging
from functools import wraps
import subprocess as sp
import uuid
import os
from os.path import join
import shutil
from ...utilities import (
    is_valid_xtb_solvent,
    XTBInvalidSolventError,
    XTBExtractor
)


logger = logging.getLogger(__name__)


class EnergyError(Exception):
    """
    Indicates a failed energy calculation.

    """

    ...


def _add_cache_use(get_energy):
    """
    Makes :meth:`~EnergyCalculator.get_energy` use the cache.

    Decorates `get_energy` so that before running it checks if the
    :class:`.Molecule` has already had its energy calculated. If so,
    and :attr:`~EnergyCalculator._use_cache` is ``True``, then the
    energy value in the cache is returned.

    Parameters
    ----------
    get_energy : :class:`function`
        A function which is to have cache use added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(get_energy)
    def inner(self, mol):
        if self._use_cache and mol in self._cache:
            logger.info(f'Using cached energy value with {mol}.')
            return self._cache[mol]
        else:
            e = get_energy(self, mol)
            if self._use_cache:
                self._cache[mol] = e
            return e

    return inner


class EnergyCalculator:
    """
    Calculates the energy of molecules.

    """

    def __init__(self, use_cache=False):
        """
        Initialize a :class:`EnergyCalculator` instance.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`get_energy` will not run twice on the
            same molecule, but will instead return the previously
            calculated value.

        """

        # Maps molecules to previously calculated energy values.
        self._cache = {}
        self._use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.get_energy = _add_cache_use(cls.get_energy)
        return super().__init_subclass__(**kwargs)

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class FormationEnergy(EnergyCalculator):
    """
    Calculates the formation energy of a molecule.

    Examples
    --------
    .. code-block:: python

        import stk

        # Make the molecules needed to calculate formation energy.
        water = stk.BuildingBlock('[H]O[H]')
        bb1 = stk.BuildingBlock('NCCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC=O', ['aldehyde'])

        chain = stk.polymer.Linear('AB', [0, 0], 6)
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=chain
        )

        # Make the energy calculator used to calculate energies.
        uff_energy = stk.UFFEnergy()

        # Make the formation energy calculator.
        formation = stk.FormationEnergy(
            energy_calculator=uff_energy,
            reactants=[bb1]*6 + [bb2]*6,
            products=[water]*polymer.bonds_made
        )

        # Calculate the formation energy.
        formation_energy = formation.get_energy(polymer)

    """

    def __init__(
        self,
        energy_calculator,
        reactants,
        products,
        use_cache=False
    ):
        """
        Initialize a :class:`FormationEnergy` instance.

        Parameters
        ----------
        energy_calculator : :class:`EnergyCalculator`
            The :class:`EnergyCalculator` used to calculate the energy
            of all reactant and product molecules.

        reactants : :class:`list` of :class:`.Molecule`
            The reactants. If there are multiples of the same reactant
            then it must appear multiple times in this :class:`list`.

        products : :class:`list` of :class:`.Molecule`
            The molecules which are produced as a result of the
            formation reaction. This :class:`list` must omit the
            :class:`.Molecule` passed to :meth:`get_energy`. If there
            are multiples of the same product, it must appear multiple
            times in this :class:`list`.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`get_energy` will not run twice on the
            same molecule, but will instead return the previously
            calculated value.

        """

        self._energy_calculator = energy_calculator
        self._reactants = reactants
        self._products = products
        super().__init__(use_cache=use_cache)

    def get_energy(self, mol):
        """
        Calculate the formation energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose formation energy is to be
            calculated.

        Returns
        -------
        :class:`float`
            The formation energy.

        """

        product_energy = sum(
            self._energy_calculator.get_energy(product)
            for product in self._products
        )
        product_energy += self._energy_calculator.get_energy(mol)

        reactant_energy = sum(
            self._energy_calculator.get_energy(reactant)
            for reactant in self._reactants
        )

        return product_energy - reactant_energy


class MMFFEnergy(EnergyCalculator):
    """
    Uses the MMFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        mmff = stk.MMFFEnergy()

        # Calculate the energy.
        energy1 = mmff.get_energy(mol1)

    """

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.MMFFGetMoleculeForceField(
            rdkit_mol,
            rdkit.MMFFGetMoleculeProperties(rdkit_mol)
        )
        return ff.CalcEnergy()


class UFFEnergy(EnergyCalculator):
    """
    Uses the UFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        uff = stk.UFFEnergy()

        # Calculate the energy.
        energy1 = uff.get_energy(mol1)

    """

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        # RingInfo needs to be initialized, else rdkit may raise an
        # error.
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.UFFGetMoleculeForceField(rdkit_mol)
        return ff.CalcEnergy()


class XTBEnergy(EnergyCalculator):
    """
    Uses GFN-xTB [1]_ to calculate energy and other properties.

    By default, :meth:`get_energy` will extract other properties of the
    :class:`.Molecule` passed to :meth:`get_energy`, which
    will be saved in the attributes of :class:`.XTBEnergy`.

    Notes
    -----
    When running :meth:`get_energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTBEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    total_energies : :class:`dict`
        :class:`dict` of the total energy of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qonly_dipole_moments : :class:`dict`
        :class:`dict` of the q only dipole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    full_dipole_moments : :class:`dict`
        :class:`dict` of the full dipole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qonly_quadrupole_moments : :class:`dict`
        :class:`dict` of the q only quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qdip_quadrupole_moments : :class:`dict`
        :class:`dict` of the q+dip quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    full_quadrupole_moments : :class:`dict`
        :class:`dict` of the full quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    total_free_energies : :class:`dict`
        :class:`dict` of the total free energy of
        each :class:`.Molecule` passed to :meth:`get_energy`.
        This is empty if :attr:`calculate_free_energy` is ``False``.

    frequencies : :class:`dict`
        :class:`dict` of the vibrational frequencies of
        each :class:`.Molecule` passed to :meth:`get_energy`.
        This is empty if `calculate_free_energy` was ``False``.

    Examples
    --------
    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCNCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear("AB", [0, 0], 3)
        )

        # Optimize the constructed molecule so that it has a
        # reasonable structure.
        optimizer = stk.OptimizerSequence(
            stk.ETKDG(),
            stk.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
        )
        optimizer.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = stk.XTBEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True
        )

        p_total_energy = xtb.get_energy(polymer)

        # Extract properties from the energy calculator for a given
        # molecule.
        homo_lumo_gap = xtb.homo_lumo_gaps[polymer]
        fermi_levels = xtb.fermi_levels[polymer]
        homo_lumo_orbitals = xtb.homo_lumo_orbitals[polymer]
        full_dipole_moments = xtb.full_dipole_moments[polymer]
        full_quadrupole_moments = xtb.full_quadrupole_moments[polymer]

        # The total energy can be extracted at any point from the
        # calculator.
        total_energy = xtb.total_energies[polymer]

    If `calculate_free_energy` is ``True``, xTB performs a
    numerical Hessian calculation and calculates the total free energy
    and vibrational frequencies of a molecule. It is recommended that a
    well optimized structure be used as input for these calculations

    .. code-block:: python

        # Optimize the constructed molecule so that it has a
        # reasonable structure.
        optimizer = stk.OptimizerSequence(
            stk.ETKDG(),
            stk.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True,
                opt_level='verytight'
            )
        )
        optimizer.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = stk.XTBEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
            calculate_free_energy=True
        )

        p_total_energy = xtb.get_energy(polymer)

        # Extract properties from the energy calculator for a given
        # molecule.
        p_total_free_energy = xtb.total_energies[polymer]
        p_frequencies = xtb.frequencies[polymer]

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html

    """
    def __init__(
        self,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        num_cores=1,
        calculate_free_energy=False,
        electronic_temperature=300,
        solvent=None,
        solvent_grid='normal',
        charge=0,
        num_unpaired_electrons=0,
        unlimited_memory=False,
        use_cache=False
    ):
        """
        Initializes a :class:`XTBEnergy` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`int`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        calculate_free_energy : :class:`bool`, optional
            Whether to calculate the total free energy and vibrational
            frequencies. Setting this to ``True`` can drastically
            increase calculation time and memory requirements.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        charge : :class:`int`, optional
            Formal molecular charge.

        num_unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints
            on the stack size. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule.

        """
        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    f'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(gfn_version, solvent):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} is invalid for ',
                    f'version {gfn_version!r}.'
                )

        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
        self._num_cores = str(num_cores)
        self._calculate_free_energy = calculate_free_energy
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory

        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.qonly_quadrupole_moments = {}
        self.qdip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        self.total_free_energies = {}
        self.frequencies = {}
        super().__init__(use_cache=use_cache)

    def _get_properties(self, mol, output_file):
        """
        Extracts properties from a GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy was calculated.

        output_file : :class: `str`
            Name of the output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Get properties from output_file.
        xtbext = XTBExtractor(output_file=output_file)

        self.total_energies[mol] = xtbext.total_energy
        self.homo_lumo_gaps[mol] = xtbext.homo_lumo_gap
        self.fermi_levels[mol] = xtbext.fermi_level
        self.homo_lumo_orbitals[mol] = xtbext.homo_lumo_occ
        self.qonly_dipole_moments[mol] = xtbext.qonly_dipole_moment
        self.full_dipole_moments[mol] = xtbext.full_dipole_moment
        self.qonly_quadrupole_moments[mol] = \
            xtbext.qonly_quadrupole_moment
        self.qdip_quadrupole_moments[mol] = \
            xtbext.qdip_quadrupole_moment
        self.full_quadrupole_moments[mol] = \
            xtbext.full_quadrupole_moment
        if self._calculate_free_energy:
            self.total_free_energies[mol] = xtbext.total_free_energy
            self.frequencies[mol] = xtbext.frequencies

    def _run_xtb(self, xyz, out_file):
        """
        Runs GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        if self._solvent is not None:
            solvent = f'--gbsa {self._solvent} {self._solvent_grid}'
        else:
            solvent = ''

        if self._calculate_free_energy:
            calc_type = '--hess'
        else:
            calc_type = ''

        cmd = (
            f'{memory} {self._xtb_path} '
            f'{xyz} --gfn {self._gfn_version} '
            f'{calc_type} --parallel {self._num_cores} '
            f'--etemp {self._electronic_temperature} '
            f'{solvent} --chrg {self._charge} '
            f'--uhf {self._num_unpaired_electrons}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        """

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        xyz = join(output_dir, 'input_structure.xyz')
        out_file = join(output_dir, 'energy.output')
        mol.write(xyz)

        try:
            os.chdir(output_dir)
            self._run_xtb(xyz=xyz, out_file=out_file)
        finally:
            os.chdir(init_dir)

        self._get_properties(
            mol=mol,
            output_file=out_file
        )
        return self.total_energies[mol]
