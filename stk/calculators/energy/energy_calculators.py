"""
Defines energy calculators.

See :mod:`.energy`.

"""

import rdkit.Chem.AllChem as rdkit
import logging
from functools import wraps
import subprocess as sp
import uuid
import os
import shutil
from ...utilities import (is_valid_xtb_solvent,
                          XTBInvalidSolventError,
                          XTBExtrators)


logger = logging.getLogger(__name__)


class EnergyError(Exception):
    """
    Indicates a failed energy calculation.

    """

    ...


def _add_cache_use(energy):
    """
    Makes :meth:`~EnergyCalculator.energy` use the cache.

    Decorates `energy` so that before running it checks if the
    :class:`.Molecule` and conformer have already had their energy
    calculated. If so, and :attr:`~EnergyCalculator.use_cache` is
    ``True``, then the energy value in the
    :attr:`~EnergyCalculator.cache` is returned.

    Parameters
    ----------
    energy : :class:`function`
        A function which is to have cache use added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(energy)
    def inner(self, mol, conformer=-1):
        key = (mol.key, conformer)
        if self.use_cache and key in self.cache:
            logger.info(
                f'Using cached energy value with '
                f'"{mol.name}" conformer {conformer}.'
            )
            return self.cache[key]
        else:
            e = energy(self, mol, conformer)
            if self.use_cache:
                self.cache[key] = e
            return e

    return inner


class EnergyCalculator:
    """
    Calculate the energy of molecules.

    Attributes
    ----------
    cache : :class:`dict`
        A :class:`dict` hodling calculated energy values of the form

        .. code-block:: python

            cache = {
                (mol1, conformer1): 12.2,
                (mol1, conformer2): 2.0,
                (mol2, conformer2): 321.12
            }

        which holds every :class:`Molecule` and conformer whose energy
        was calculated by the :class:`EnergyCalculator`. Here ``mol1``
        and ``mol2`` are :class:`.Molecule` objects and ``conformer1``
        and ``conformer2`` are of type :class:`int`, and are the ids
        of the conformers passed to :meth:`energy`.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer, but will instead return the previously
        calculated value.

    """

    def __init__(self, use_cache=False):
        """
        Initializes a :class:`EnergyCalculator` instance.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer, but will instead return the
            previously calculated value.

        """

        self.cache = {}
        self.use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.energy = _add_cache_use(cls.energy)
        return super().__init_subclass__(**kwargs)

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        raise NotImplementedError()


class FormationEnergy(EnergyCalculator):
    """
    Calculates the formation energy of a molecule.

    Attributes
    ----------
    energy_calculator : :class:`EnergyCalculator`
        The :class:`EnergyCalculator` used to calculate the energy of
        all reactant and product molecules.

    reactants : :class:`list` of :class:`.Molecule`
        The reactants. If there are multiples of the same reactant
        then it must appear multiple times in this :class:`list`.

    products : :class:`list` of :class:`.Molecule`
        The molecules which are produced as a result of the formation
        reaction. This :class:`list` must omit the :class:`.Molecule`
        passed to :meth:`energy`. If there are multiples of the same
        product, it must appear multiple times in this :class:`list`.

    reactant_conformers : :class:`list` of :class:`int`
        The conformer ids of the molecules in :attr:`reactants` to use.

    product_conformers : :class:`list` of :class:`int`
        The conformer ids of the molecules in :attr:`products` to use.

    Examples
    --------
    .. code-block:: python

        # Make the molecules needed to calculate formation energy.
        water = StructUnit.smiles_init('[H]O[H]')
        bb1 = StructUnit2(...)
        bb2 = StructUnit3(...)
        cage = Cage([bb1, bb2], FourPlusSix())

        # Make the energy calculator used to calculate energies.
        uff_energy = UFFEnergy()

        # Make the formation energy calculator.
        formation = FormationEnergy(
            energy_calculator=uff_energy,
            reactants=[bb1]*6 + [bb2]*4,
            products=[water]*cage.topology.bonds_made
        )

        # Calculate the formation energy.
        formation_energy = formation.energy(cage)

    """

    def __init__(self,
                 energy_calculator,
                 reactants,
                 products,
                 reactant_conformers=None,
                 product_conformers=None,
                 use_cache=False):
        """
        Initializes a :class:`FormationEnergy` instance.

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
            :class:`.Molecule` passed to :meth:`energy`. If there are
            multiples of the same product, it must appear multiple
            times in this :class:`list`.

        reactant_conformers : :class:`list` of :class:`int`, optional
            The conformer ids of the molecules in :attr:`reactants` to
            use.

        product_conformers : :class:`list` of :class:`int`, optional
            The conformer ids of the molecules in :attr:`products` to
            use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer, but will instead return the
            previously calculated value.

        """

        if reactant_conformers is None:
            reactant_conformers = [-1 for i in range(len(reactants))]

        if product_conformers is None:
            product_conformers = [-1 for i in range(len(products))]

        self.energy_calculator = energy_calculator
        self.reactants = reactants
        self.products = products
        self.reactant_conformers = reactant_conformers
        self.product_conformers = product_conformers
        super().__init__(use_cache=use_cache)

    def energy(self, mol, conformer=-1):
        """
        Calculates the formation energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose formation energy is to be
            calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The formation energy.

        """

        products = zip(self.products, self.product_conformers)
        product_energy = sum(
            self.energy_calculator.energy(product, conformer)
            for product, conformer in products
        )
        product_energy += self.energy_calculator.energy(mol, conformer)

        reactants = zip(self.reactants, self.reactant_conformers)
        reactant_energy = sum(
            self.energy_calculator.energy(reactant, conformer)
            for reactant, conformer in reactants
        )

        return product_energy - reactant_energy


class MMFFEnergy(EnergyCalculator):
    """
    Uses the MMFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        # Create a molecules whose energy we want to know.
        mol1 = StructUnit.smiles_init('CCCNCCCN')
        mol2 = Polymer(...)
        mol3 = Cage(...)

        # Create the energy calculator.
        mmff = MMFFEnergy()

        # Calculate the energies.
        energy1 = mmff.energy(mol1)
        energy2 = mmff.energy(mol2)
        energy3 = mmff.energy(mol3)

    """

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit.GetSSSR(mol.mol)
        mol.mol.UpdatePropertyCache()
        ff = rdkit.MMFFGetMoleculeForceField(
            mol.mol,
            rdkit.MMFFGetMoleculeProperties(mol.mol),
            confId=conformer
        )

        return ff.CalcEnergy()


class UFFEnergy(EnergyCalculator):
    """
    Uses the UFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        # Create a molecules whose energy we want to know.
        mol1 = StructUnit.smiles_init('CCCNCCCN')
        mol2 = Polymer(...)
        mol3 = Cage(...)

        # Create the energy calculator.
        uff = UFFEnergy()

        # Calculate the energies.
        energy1 = uff.energy(mol1)
        energy2 = uff.energy(mol2)
        energy3 = uff.energy(mol3)

    """

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optinal
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        mol.mol.UpdatePropertyCache()
        # RingInfo needs to be initialized, else rdkit may raise an
        # error.
        rdkit.GetSSSR(mol.mol)
        ff = rdkit.UFFGetMoleculeForceField(mol.mol, confId=conformer)
        return ff.CalcEnergy()


class XTBEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate energy and other properties.

    By default, :meth:`energy` will extract other properties of the
    :class:`.Molecule` and conformer passed to :meth:`energy`, which
    will be saved in attributes of :class:`.XTBEnergy`.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTBEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`int`
        Parameterization of GFN to use in xTB.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/basics.html.

    output_dir : :class:`str`
        The name of the directory into which files generated during
        the calculation are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores xTB should use.

    electronic_temperature : :class:`int`
        Electronic temperature to use (in K).

    solvent : :class:`str`
        Solvent to use in GBSA implicit solvation method.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    solvent_grid : :class:`str`
        Grid level to use in SASA calculations for GBSA implicit
        solvent.
        Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
        or ``'extreme'``.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    charge : :class:`int`
        Formal molecular charge.

    num_unpaired_electrons : :class:`int`
        Number of unpaired electrons.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    unlimited_memory : :class: `bool`
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    total_energies : :class:`dict`
        :class:`dict` of the total energy of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of :class:`.Molecule` and
        conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qonly_dipole_moments : :class:`dict`
        :class:`dict` of the `q only` dipole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    full_dipole_moments : :class:`dict`
        :class:`dict` of the `full` dipole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qonly_quadrupole_moments : :class:`dict`
        :class:`dict` of the `q only` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qdip_quadrupole_moments : :class:`dict`
        :class:`dict` of the `q+dip` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    full_quadrupole_moments : :class:`dict`
        :class:`dict` of the `full` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    Examples
    --------

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBEnergy` should usually be used after optimization with
    some other method. This is because xTB only uses xyz coordinates
    as input and so will not recognize the long bonds created during
    assembly. An optimizer which can minimize these bonds should be
    used before :class:`XTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        # Pre optimize molecule with UFF.
        uff = UFF()
        uff.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = XTBEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
        )

        # Runs the calculation and returns energy.
        p_totalenergy = xtb.energy(polymer, conformer)

        # Extracts properties from energy calculator for a given
        # molecule and conformer.
        key = (polymer, conformer)
        p_homo_lumo_gap = xtb.homo_lumo_gaps[key]
        p_fermi_levels = xtb.fermi_levels[key]
        p_homo_lumo_orbitals = xtb.homo_lumo_orbitals[key]
        p_full_dipole_moments = xtb.full_dipole_moments[key]
        p_full_quadrupole_moments = xtb.full_quadrupole_moments[key]

        # The total energy can be extracted at any point from the
        # calculator.
        polymer_totalenergy = xtb.total_energies[key]

    See Also
    --------
    #. https://xtb-docs.readthedocs.io/en/latest/setup.html

    """
    def __init__(self,
                 xtb_path,
                 gfn_version=2,
                 output_dir=None,
                 num_cores=1,
                 electronic_temperature=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=0,
                 num_unpaired_electrons=0,
                 use_cache=False,
                 unlimited_memory=False):
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

        electronic_temperature : :class:`int`, optional
            Electronic temperature to use (in K).

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

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints
            on the stacksize. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        """
        if solvent is not None:
            solvent = solvent.lower()
            if not is_valid_xtb_solvent(gfn_version, solvent):
                if gfn_version == '0':
                    raise XTBInvalidSolventError(
                        f'No solvent valid for version',
                        f' {gfn_version!r}.'
                    )
                else:
                    raise XTBInvalidSolventError(
                        f'Solvent {solvent!r} is invalid for ',
                        f'version {gfn_version!r}.'
                    )

        self.xtb_path = xtb_path
        self.gfn_version = str(gfn_version)
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.electronic_temperature = str(electronic_temperature)
        self.solvent = solvent
        self.solvent_grid = solvent_grid
        self.charge = str(charge)
        self.num_unpaired_electrons = str(num_unpaired_electrons)
        self.unlimited_memory = unlimited_memory

        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.qonly_quadrupole_moments = {}
        self.qdip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        super().__init__(use_cache=use_cache)

    def _get_properties(self, mol, conformer, output_file):
        """
        Extracts properties from a GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`
            The conformer of `mol` to use.

        output_file : :class: `str`
            Name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Get properties from output_file.
        xtbext = XTBExtrators(output_file=output_file)
        key = (mol, conformer)
        self.total_energies[key] = xtbext.total_energy()
        self.homo_lumo_gaps[key] = xtbext.homo_lumo_gap()
        self.fermi_levels[key] = xtbext.fermi_level()
        self.homo_lumo_orbitals[key] = xtbext.homo_lumo_occ()
        self.qonly_dipole_moments[key] = xtbext.qonly_dipole_moment()
        self.full_dipole_moments[key] = xtbext.full_dipole_moment()
        self.qonly_quadrupole_moments[key] = \
            xtbext.qonly_quadrupole_moment()
        self.qdip_quadrupole_moments[key] = \
            xtbext.qdip_quadrupole_moment()
        self.full_quadrupole_moments[key] = \
            xtbext.full_quadrupole_moment()

    def _write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN-xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`
            The conformer of `mol` to use.

        Returns
        -------
        :class:`str`
            Name of output file with xTB results.

        """
        xyz = 'input_structure.xyz'
        out_file = 'energy.output'
        mol.write(xyz, conformer=conformer)
        # Modify the memory limit.
        if self.unlimited_memory:
            # Uses the shell if unlimited_memory is True to be run
            # multiple commpands in one subprocess.
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        if self.solvent is not None:
            if self.solvent_grid == 'normal':
                solvent_grid = ''
            else:
                solvent_grid = ''
            solvent = f'--gbsa {self.solvent} {solvent_grid}'
        else:
            solvent = ''

        cmd = (
            f'{memory} {self.xtb_path} {xyz} --gfn {self.gfn_version} '
            f'--parallel {self.num_cores} '
            f'--etemp {self.electronic_temperature} '
            f'{solvent} --chrg {self.charge} '
            f'--uhf {self.num_unpaired_electrons}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                shell=self.unlimited_memory
            )

        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optional
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        try:
            os.chdir(output_dir)
            out_file = self._write_and_run_command(
                mol=mol,
                conformer=conformer
            )
            self._get_properties(
                mol=mol,
                conformer=conformer,
                output_file=out_file
            )
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]


class XTBFreeEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate free energy and other properties.

    By default, :meth:`energy` will extract other properties of the
    :class:`.Molecule` and conformer passed to :meth:`energy`, which
    will be saved in attributes of :class:`.XTBFreeEnergy`.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTBFreeEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`int`
        Parameterization of GFN to use in xTB.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/basics.html.

    output_dir : :class:`str`
        The name of the directory into which files generated during
        the calculation are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores xTB should use.

    electronic_temperature : :class:`int`
        Electronic temperature to use (in K).

    solvent : :class:`str`
        Solvent to use in GBSA implicit solvation method.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    solvent_grid : :class:`str`
        Grid level to use in SASA calculations for GBSA implicit
        solvent.
        Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
        or ``'extreme'``.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    charge : :class:`int`
        Formal molecular charge.

    num_unpaired_electrons : :class:`int`
        Number of unpaired electrons.

    use_cache : :class:`bool`
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    unlimited_memory : :class: `bool`
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be
        ``True``, however this may raise issues on clusters.

    total_energies : :class:`dict`
        :class:`dict` of the total energy of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qonly_dipole_moments : :class:`dict`
        :class:`dict` of the `q only` dipole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    full_dipole_moments : :class:`dict`
        :class:`dict` of the `full` dipole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qonly_quadrupole_moments : :class:`dict`
        :class:`dict` of the `q only` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    qdip_quadrupole_moments : :class:`dict`
        :class:`dict` of the `q+dip` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    full_quadrupole_moments : :class:`dict`
        :class:`dict` of the `full` quadrupole momentent of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    total_free_energies : :class:`dict`
        :class:`dict` of the total free energy of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.

    frequencies : :class:`dict`
        :class:`dict` of the vibrational frequencies of
        :class:`.Molecule` and conformer passed to :meth:`energy`.
        A :class:`dict` mapping each :class:`.Molecule` and Key has the
        form ``(mol, conformer)``.


    Examples
    --------

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBFreeEnergy` should usually be used after optimization
    with some other method. This is because xTB only uses xyz
    coordinates as input and so will not recognize the long bonds
    created during assembly. An optimizer which can minimize these
    bonds should be used before :class:`XTBFreeEnergy`.

    :class:`XTBFreeEnergy` performs a numerical Hessian calculation
    and calculates the total free energy and vibrational frequencies
    of a molecule. It is recommended that a well optimized structure
    be used as input for these calculations. For example, very tight
    criteria in the :class:`XTB` optimizer will provide a reasonable
    structure for this analysis.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        # Pre optimize molecule with UFF.
        uff = UFF()

        # Optimization with XTB using very tight convergence
        # criteria.
        xtb_opt = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
            opt_level='verytight',
        )

        # Define optimizer sequence and run optimzations.
        opt = OptimizerSequence(uff, xtb_opt)
        opt.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = XTBFreeEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
        )

        # Runs the calculation and returns energy.
        p_totalenergy = xtb.energy(polymer, conformer)

        # Extracts properties from energy calculator for a given
        # molecule and conformer.
        key = (polymer, conformer)
        # Extract the total free energy and vibrational frequencies.
        p_total_free_energy = xtb.total_energies[key]
        p_frequencies = xtb.frequencies[key]

        # All the properties available in :class:`XTBEnergy` are also
        # available.

    See Also
    --------
    #. https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(self,
                 xtb_path,
                 gfn_version=2,
                 output_dir=None,
                 num_cores=1,
                 electronic_temperature=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=0,
                 num_unpaired_electrons=0,
                 use_cache=False,
                 unlimited_memory=False):
        """
        Initializes a :class:`XTBFreeEnergy` instance.

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

        electronic_temperature : :class:`int`, optional
            Electronic temperature to use (in K).

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

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

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints
            on the stacksize. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        """
        if solvent is not None:
            solvent = solvent.lower()
            if not is_valid_xtb_solvent(gfn_version, solvent):
                if gfn_version == '0':
                    raise XTBInvalidSolventError(
                        f'No solvent valid for version',
                        f' {gfn_version!r}.'
                    )
                else:
                    raise XTBInvalidSolventError(
                        f'Solvent {solvent!r} is invalid for ',
                        f'version {gfn_version!r}.'
                    )

        self.xtb_path = xtb_path
        self.gfn_version = str(gfn_version)
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.electronic_temperature = str(electronic_temperature)
        self.solvent = solvent
        self.solvent_grid = solvent_grid
        self.charge = str(charge)
        self.num_unpaired_electrons = str(num_unpaired_electrons)
        self.unlimited_memory = unlimited_memory

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

    def _get_properties(self, mol, conformer, output_file):
        """
        Extracts properties from a GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`
            The conformer of `mol` to use.

        output_file : :class: `str`
            Name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Get properties from output_file.
        xtbext = XTBExtrators(output_file=output_file)
        key = (mol, conformer)
        self.total_energies[key] = xtbext.total_energy()
        self.homo_lumo_gaps[key] = xtbext.homo_lumo_gap()
        self.fermi_levels[key] = xtbext.fermi_level()
        self.homo_lumo_orbitals[key] = xtbext.homo_lumo_occ()
        self.qonly_dipole_moments[key] = xtbext.qonly_dipole_moment()
        self.full_dipole_moments[key] = xtbext.full_dipole_moment()
        self.qonly_quadrupole_moments[key] = \
            xtbext.qonly_quadrupole_moment()
        self.qdip_quadrupole_moments[key] = \
            xtbext.qdip_quadrupole_moment()
        self.full_quadrupole_moments[key] = \
            xtbext.full_quadrupole_moment()
        self.total_free_energies[key] = xtbext.total_free_energy()
        self.frequencies[key] = xtbext.frequencies()

    def _write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN-xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`
            The conformer of `mol` to use.

        Returns
        -------
        :class:`str`
            Name of output file with xTB results.

        """
        xyz = 'input_structure.xyz'
        out_file = 'free_energy.output'
        mol.write(xyz, conformer=conformer)
        # Modify the memory limit.
        if self.unlimited_memory:
            # Uses the shell if unlimited_memory is True to be run
            # multiple commpands in one subprocess.
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        if self.solvent is not None:
            if self.solvent_grid == 'normal':
                solvent_grid = ''
            else:
                solvent_grid = ''
            solvent = f'--gbsa {self.solvent} {solvent_grid}'
        else:
            solvent = ''

        cmd = (
            f'{memory} {self.xtb_path} {xyz} --gfn {self.gfn_version} '
            f'--hess --parallel {self.num_cores} '
            f'--etemp {self.electronic_temperature} '
            f'{solvent} --chrg {self.charge} '
            f'--uhf {self.num_unpaired_electrons}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                shell=self.unlimited_memory
            )

        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        conformer : :class:`int`, optional
            The conformer of `mol` to use.

        Returns
        -------
        :class:`float`
            The energy.

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        try:
            os.chdir(output_dir)
            out_file = self._write_and_run_command(
                mol=mol,
                conformer=conformer
            )
            self._get_properties(
                mol=mol,
                conformer=conformer,
                output_file=out_file
            )
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]
