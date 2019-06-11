"""
Defines energy calculators.

See :mod:`.energy`.

"""

import rdkit.Chem.AllChem as rdkit
import logging
import re
from functools import wraps
import subprocess as sp
import uuid
import os
import shutil
from ...utilities import valid_XTB_solvent


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
            molecule and conformer, but will instead return the previously
            calculated value.

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


class XTBEnergyHessianFailedError(Exception):
    ...


class XTBEnergyNegativeFreqError(Exception):
    ...


class XTBEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate energies and other properties of molecules.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.EnergyCalculator` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Documentation for GFN2-xTB available:
    https://xtb-docs.readthedocs.io/en/latest/setup.html

    Attributes
    ----------
    gfnxtb_path : :class:`str`
        The path to the GFN-xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN-xTB to use.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores for GFN-xTB to use. Requires appropriate setup
        of GFN-xTB by user.

    use_cache : :class:`bool`, optional
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    mem_ulimit : :class: `bool`, optional
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be ``True``,
        however this may raise issues on clusters.

    etemp : :class:`int`, optional
        Electronic temperature to use (in K). Defaults to 300K.

    solvent : :class:`str`, optional
        Solvent to use in GBSA implicit solvation method.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`, optional
        Grid level to use in SASA calculations for GBSA implicit solvent.
        Options:
            normal, tight, verytight, extreme
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    multiplicity : :class:`str`, optional
        Number of unpaired electrons.

    charge : :class:`str`, optional
        Formal molecular charge. `-` should be used to indicate sign.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = XTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBEnergy` should usually be used after optimization with some
    other method. This is because GFN-xTB only uses xyz coordinates as input
    and so will not recognize the long bonds created during assembly.
    An optimizer which can minimize these bonds should be used before
    :class:`XTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        gfnxtb = XTBEnergy(gfnxtb_path='/opt/gfnxtb/xtb',
                              mem_ulimit=True)
        gfnxtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBEnergy`. For example vibrational frequencies,
    the HOMO-LUMO gap and thermodynamic properties (such as the Total Free
    Energy) can be calculated after an optimization with very tight contraints
    with an implicit solvent (THF). Very tight criteria are required to ensure
    that no negative vibrational frequencies are present.

    .. code-block:: python

        gfnxtb = OptimizerSequence(
            UFF(),
            XTB(gfnxtb_path='/opt/gfnxtb/xtb',
                   mem_ulimit=True,
                   opt_level='verytight',
                   solvent='THF')
        )
        gfnxtb.optimize(polymer)

        gfnxtbenergy = XTBEnergy(gfnxtb_path='/opt/gfnxtb/xtb',
                                    mem_ulimit=True,
                                    free=True,
                                    solvent='THF')
        polymer_properties = gfnxtbenergy.energy(polymer)
        polymer_free_energy = polymer_properties['totalfreeenergy']
        polymer_freq = polymer_properties['frequencies']
        polymer_gap = polymer_properties['HLGap']
    """
    def __init__(self,
                 gfnxtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 multiplicity=None,
                 use_cache=False,
                 mem_ulimit=False):
        """
        Initializes a :class:`XTBEnergy` instance.

        Parameters
        ----------
        gfnxtb_path : :class:`str`
            The path to the GFN-xTB or GFN2-xTB executable.

        gfn_version : :class:`str`
            Parameterization of GFN-xTB to use.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`
            The number of cores for GFN-xTB to use. Requires appropriate setup
            of GFN-xTB by user.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`energy` will not run twice on the same
            molecule and conformer.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`energy` will be run without constraints on
            the stacksize. If memory issues are encountered, this should be ``True``,
            however this may raise issues on clusters.

        etemp : :class:`int`, optional
            Electronic temperature to use (in K). Defaults to 300K.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit solvent.
            Options:
                normal, tight, verytight, extreme
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        multiplicity : :class:`str`, optional
            Number of unpaired electrons.

        charge : :class:`str`, optional
            Formal molecular charge. `-` should be used to indicate sign.

        """
        self.gfnxtb_path = gfnxtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.num_cores = str(num_cores)
        self.etemp = str(etemp)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_XTB_solvent(gfn_version=self.gfn_version,
                                 solvent=self.solvent)
        self.solvent_grid = solvent_grid
        self.charge = charge
        self.multiplicity = multiplicity
        self.mem_ulimit = mem_ulimit
        # properties
        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.Qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.Qonly_quadrupole_moments = {}
        self.QDip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        super().__init__(use_cache=use_cache)

    def __ext_total_energy(self, output_string):
        """
        Extracts total energy (a.u.) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "          | TOTAL ENERGY              -76.260405590154 Eh   |"

        Returns
        -------
        :class:`float`
            Total energy in a.u.
        """
        value = None

        # regex for numbers
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        for line in reversed(output_string):
            if '          | TOTAL ENERGY  ' in line:
                value = nums.search(line.rstrip()).group(0)
                break

        return float(value)

    def __ext_homo_lumo_gap(self, output_string):
        """
        Extracts total energy (eV) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "          | HOMO-LUMO GAP               2.336339660160 eV   |"

        Returns
        -------
        :class:`float`
            Homo-Lumo gap in eV.
        """
        value = None

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        for line in reversed(output_string):
            if '          | HOMO-LUMO GAP   ' in line:
                value = nums.search(line.rstrip()).group(0)
                break

        return float(value)

    def __ext_fermi_level(self, output_string):
        """
        Extracts Fermi-Level energy (eV) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "             Fermi-level           -0.3159871 Eh           -8.5984 eV"

        Returns
        -------
        :class:`float`
            Fermi-level in eV.
        """
        value = None

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        for line in reversed(output_string):
            if '             Fermi-level        ' in line:
                part2 = line.split('Eh')
                value = nums.search(part2[1].rstrip()).group(0)
                break

        return float(value)

    def __ext_Qonly_dipole_emom(self, output_string):
        """
        Extracts `q only` dipole moment vector (Debye) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        " q only:       -0.033      -0.081      -0.815"

        Returns
        -------
        :class:`list` of :class:`float`
            Components of dipole moment in a list of length 3.
            [`x`, `y`, `z`]
        """
        value = None

        sample_set = []
        for i, line in enumerate(output_string):
            if 'molecular dipole:' in line:
                sample_set = output_string[i+2].rstrip()

        # get values from line
        if 'q only:' in sample_set:
            x, y, z = [i for i in sample_set.split(':')[1].split(' ')
                       if i != '']

        value = [float(x), float(y), float(z)]

        return value

    def __ext_full_dipole_mom(self, output_string):
        """
        Extracts `full` dipole moment vector (Debye) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "   full:       -0.684       0.122      -1.071       3.245"

        Returns
        -------
        :class:`list` of :class:`float`
            Components of dipole moment and total magnitude in a list of length 4.
            [`x`, `y`, `z`, `tot (Debye)`]
        """
        value = None

        sample_set = []
        for i, line in enumerate(output_string):
            if 'molecular dipole:' in line:
                sample_set = output_string[i+3].rstrip()

        # get values from line
        if 'full:' in sample_set:
            x, y, z, m = [i for i in sample_set.split(':')[1].split(' ')
                          if i != '']

        value = [float(x), float(y), float(z), float(m)]

        return value

    def __ext_Qonly_quadrupole_mom(self, output_string):
        """
        Extracts `q only` traceless quadrupole moment vector (Debye) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        " q only:        7.152      10.952       3.364      15.349       2.074     -10.515"

        Returns
        -------
        :class:`list` of :class:`float`
            Components of quadrupole moment in a list of length 6.
            [`xx`, `xy`, `xy`, `xz`, `yz`, `zz`]
        """
        value = None

        sample_set = []
        for i, line in enumerate(output_string):
            if 'molecular quadrupole (traceless):' in line:
                sample_set = output_string[i+2].rstrip()

        # get values from line
        if 'q only:' in sample_set:
            xx, xy, yy, xz, yz, zz = [i for i in sample_set.split(':')[1].split(' ')
                                      if i != '']

        value = [float(xx), float(xy), float(yy), float(xz), float(yz), float(zz)]

        return value

    def __ext_QDip_quadrupole_mom(self, output_string):
        """
        Extracts `q+dip` traceless quadrupole moment vector (Debye) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "  q+dip:       -6.239      21.552      16.601      12.864       2.504     -10.362"

        Returns
        -------
        :class:`list` of :class:`float`
            Components of quadrupole moment in a list of length 6.
            [`xx`, `xy`, `xy`, `xz`, `yz`, `zz`]
        """
        value = None

        sample_set = []
        for i, line in enumerate(output_string):
            if 'molecular quadrupole (traceless):' in line:
                sample_set = output_string[i+3].rstrip()

        # get values from line
        if 'q+dip:' in sample_set:
            xx, xy, yy, xz, yz, zz = [i for i in sample_set.split(':')[1].split(' ')
                                      if i != '']

        value = [float(xx), float(xy), float(yy), float(xz), float(yz), float(zz)]

        return value

    def __ext_full_quadrupole_mom(self, output_string):
        """
        Extracts `full` traceless quadrupole moment vector (Debye) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "   full:       -6.662      22.015      16.959      12.710       3.119     -10.297"

        Returns
        -------
        :class:`list` of :class:`float`
            Components of quadrupole moment in a list of length 6.
            [`xx`, `xy`, `xy`, `xz`, `yz`, `zz`]
        """
        value = None

        sample_set = []
        for i, line in enumerate(output_string):
            if 'molecular quadrupole (traceless):' in line:
                sample_set = output_string[i+4].rstrip()

        # get values from line
        if 'full:' in sample_set:
            xx, xy, yy, xz, yz, zz = [i for i in sample_set.split(':')[1].split(' ')
                                      if i != '']

        value = [float(xx), float(xy), float(yy), float(xz), float(yz), float(zz)]

        return value

    def __ext_homo_lumo_occ(self, output_string):
        """
        Extracts Orbital Energies and Occupations (eV) of the HOMO and LUMO from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "        70        2.0000           -0.3514143              -9.5625 (HOMO)"
        "        71                         -0.2712405              -7.3808 (LUMO)"

        Returns
        -------
        :class:`dict`
            Dictionary of (#, occupation, Energy (eV)) of HOMO and LUMO orbital
        """
        value = None

        for line in reversed(output_string):
            if '(HOMO)' in line:
                split_line = [i for i in line.rstrip().split(' ') if i != '']
                # line is: Number, occupation, energy (Ha), energy (ev), label
                # keep: Number, occupation, energy (eV)
                homo_val = [int(split_line[0]), float(split_line[1]), float(split_line[3])]
            if '(LUMO)' in line:
                split_line = [i for i in line.rstrip().split(' ') if i != '']
                # line is: Number, energy (Ha), energy (ev), label
                # keep: Number, energy (eV)
                lumo_val = [int(split_line[0]), float(0), float(split_line[2])]

        value = {'HOMO': homo_val,
                 'LUMO': lumo_val}

        return value

    def __get_properties(self, mol, conformer, output_file):
        """
        Extracts desired properties from GFN-xTB single point energy calculation.

        """
        # get output file in string
        output_string = open(output_file, 'r').readlines()

        # get properties from output string
        self.total_energies[(mol, conformer)] = self.__ext_total_energy(output_string)
        self.homo_lumo_gaps[(mol, conformer)] = self.__ext_homo_lumo_gap(output_string)
        self.fermi_levels[(mol, conformer)] = self.__ext_fermi_level(output_string)
        self.homo_lumo_orbitals[(mol, conformer)] = self.__ext_homo_lumo_occ(output_string)
        self.Qonly_dipole_moments[(mol, conformer)] = self.__ext_Qonly_dipole_emom(output_string)
        self.full_dipole_moments[(mol, conformer)] = self.__ext_full_dipole_mom(output_string)
        self.Qonly_quadrupole_moments[(mol, conformer)] = self.__ext_Qonly_quadrupole_mom(output_string)
        self.QDip_quadrupole_moments[(mol, conformer)] = self.__ext_QDip_quadrupole_mom(output_string)
        self.full_quadrupole_moments[(mol, conformer)] = self.__ext_full_quadrupole_mom(output_string)

    def __write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN.

        Returns
        -------
        out_file : :class:`str`
            Name of output file with GFN-xTB results.
        """
        xyz = 'input_structure.xyz'
        out_file = 'energy.output'
        mol.write(xyz, conformer=conformer)
        # modify memory limit
        if self.mem_ulimit:
            cmd = ['ulimit -s unlimited ;']
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.gfnxtb_path)
        cmd.append(xyz)
        # set GFN Parameterization
        if self.gfn_version != '2':
            cmd.append('--gfn')
            cmd.append(self.gfn_version)
        # set number of cores
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # add eletronic temp term
        if self.etemp != '300':
            cmd.append('--etemp')
            cmd.append(self.etemp)
        # write solvent section of cmd
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # write charge section of cmd
        if self.charge is not None:
            cmd.append('--chrg')
            cmd.append(self.charge)
        # write multiplicity section of cmd
        if self.multiplicity is not None:
            cmd.append('--uhf')
            cmd.append(self.multiplicity)
        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # uses the shell if mem_ulimit = True and waits until
        # subproces is complete. This is required to run the mem_ulimit_cmd
        # and GFN calculation in one command, which is then closed, which
        # minimizes the risk of unrestricting the memory limits.
        sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                 shell=shell)
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of molecule `mol` using GFN-xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose energy should be claculated.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        self.total_energies[(mol, conformer)] : :class:`float`
            Total energy of system.
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
            out_file = self.__write_and_run_command(mol=mol,
                                                    conformer=conformer)
            self.__get_properties(mol=mol,
                                  conformer=conformer,
                                  output_file=out_file)
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]


class XTBFreeEnergy(XTBEnergy):
    """
    Uses GFN-xTB to calculate free energies, vibrational frequencies and
    other properties of molecules.

    Notes
    -----
    When running :meth:`energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.EnergyCalculator` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Documentation for GFN2-xTB available:
    https://xtb-docs.readthedocs.io/en/latest/setup.html

    Attributes
    ----------
    gfnxtb_path : :class:`str`
        The path to the GFN-xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN-xTB to use.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores for GFN-xTB to use. Requires appropriate setup
        of GFN-xTB by user.

    use_cache : :class:`bool`, optional
        If ``True`` :meth:`energy` will not run twice on the same
        molecule and conformer.

    mem_ulimit : :class: `bool`, optional
        If ``True`` :meth:`energy` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be ``True``,
        however this may raise issues on clusters.

    etemp : :class:`int`, optional
        Electronic temperature to use (in K). Defaults to 300K.

    solvent : :class:`str`, optional
        Solvent to use in GBSA implicit solvation method.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`, optional
        Grid level to use in SASA calculations for GBSA implicit solvent.
        Options:
            normal, tight, verytight, extreme
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    charge : :class:`str`, optional
        Formal molecular charge. `-` should be used to indicate sign.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = XTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTBEnergy` should usually be used after optimization with some
    other method. This is because GFN-xTB only uses xyz coordinates as input
    and so will not recognize the long bonds created during assembly.
    An optimizer which can minimize these bonds should be used before
    :class:`XTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        gfnxtb = XTBEnergy(gfnxtb_path='/opt/gfnxtb/xtb',
                              mem_ulimit=True)
        gfnxtb.energy(polymer)

    Energies and other properties of optimized structures can be
    extracted using :class:`XTBEnergy`. For example vibrational frequencies,
    the HOMO-LUMO gap and thermodynamic properties (such as the Total Free
    Energy) can be calculated after an optimization with very tight contraints
    with an implicit solvent (THF). Very tight criteria are required to ensure
    that no negative vibrational frequencies are present.

    .. code-block:: python

        gfnxtb = OptimizerSequence(
            UFF(),
            XTB(gfnxtb_path='/opt/gfnxtb/xtb',
                   mem_ulimit=True,
                   opt_level='verytight',
                   solvent='THF')
        )
        gfnxtb.optimize(polymer)

        gfnxtbenergy = XTBEnergy(gfnxtb_path='/opt/gfnxtb/xtb',
                                    mem_ulimit=True,
                                    free=True,
                                    solvent='THF')
        polymer_properties = gfnxtbenergy.energy(polymer)
        polymer_free_energy = polymer_properties['totalfreeenergy']
        polymer_freq = polymer_properties['frequencies']
        polymer_gap = polymer_properties['HLGap']
    """

    def __init__(self,
                 gfnxtb_path,
                 gfn_version='2',
                 output_dir=None,
                 num_cores=1,
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 use_cache=False,
                 mem_ulimit=False):
        XTBEnergy.__init__(self,
                           gfnxtb_path,
                           gfn_version,
                           output_dir,
                           num_cores,
                           etemp,
                           solvent,
                           solvent_grid,
                           charge,
                           use_cache,
                           mem_ulimit)
        self.total_free_energies = {}
        self.frequencies = {}

    def __ext_total_free_energy(self, output_string):
        """
        Extracts total free energy (a.u.) from GFN-xTB output at T=298.15K.

        Formatting based on latest version of GFN-xTB (190418)
        Example line:
        "          | TOTAL FREE ENERGY         -75.832501154309 Eh   |"

        Returns
        -------
        :class:`float`
            Total free energy in a.u.
        """
        # check that hessian was performed on geometry optimized structure
        # raise error if not
        check_hessian = True
        for line in reversed(output_string):
            if '#WARNING! Hessian on incompletely optimized geometry!' in line:
                check_hessian = False
                break
        if check_hessian is False:
            raise XTBEnergyHessianFailedError(
                f'Hessian calculation performed on unoptimized structure.'
            )
        value = None

        # regex for numbers
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        for line in reversed(output_string):
            if '          | TOTAL FREE ENERGY  ' in line:
                value = nums.search(line.rstrip()).group(0)
                break

        return float(value)

    def __ext_frequencies(self, output_string):
        """
        Extracts projected vibrational frequencies (cm-1) from GFN-xTB output.

        Formatting based on latest version of GFN-xTB (190418).
        Example line:
        "eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00"

        Returns
        -------
        :class:`list`
            List of all vibrational frequencies as :class:`float`
        """
        value = None

        # use a switch to make sure we are extracting values after the
        # final property readout
        switch = False

        frequencies = []
        for i, line in enumerate(output_string):
            if '|               Frequency Printout                |' in line:
                # turn on reading as final frequency printout has begun
                switch = True
            if ' reduced masses (amu)' in line:
                # turn off reading as frequency section is done
                switch = False
            if 'eigval :' in line and switch is True:
                split_line = [i for i in line.rstrip().split(':')[1].split(' ')
                              if i != '']
                for freq in split_line:
                    frequencies.append(freq)

        value = [float(i) for i in frequencies]

        if min(value) < 0:
            raise XTBEnergyNegativeFreqError(
                'Negative frequency encountered.'
                'Structures should be optimized prior to free energy calculation.'
            )

        return value

    def __get_properties(self, mol, conformer, output_file):
        """
        Extracts desired properties from GFN-xTB single point energy calculation.

        """
        # get output file in string
        output_string = open(output_file, 'r').readlines()

        # get properties from output string
        self.total_energies[(mol, conformer)] = self.__ext_total_energy(output_string)
        self.total_free_energies[(mol, conformer)] = self.__ext_total_free_energy(output_string)
        self.frequencies[(mol, conformer)] = self.__ext_frequencies(output_string)
        self.homo_lumo_gaps[(mol, conformer)] = self.__ext_homo_lumo_gap(output_string)
        self.fermi_levels[(mol, conformer)] = self.__ext_fermi_level(output_string)
        self.homo_lumo_orbitals[(mol, conformer)] = self.__ext_homo_lumo_occ(output_string)
        self.Qonly_dipole_moments[(mol, conformer)] = self.__ext_Qonly_dipole_emom(output_string)
        self.full_dipole_moments[(mol, conformer)] = self.__ext_full_dipole_mom(output_string)
        self.Qonly_quadrupole_moments[(mol, conformer)] = self.__ext_Qonly_quadrupole_mom(output_string)
        self.QDip_quadrupole_moments[(mol, conformer)] = self.__ext_QDip_quadrupole_mom(output_string)
        self.full_quadrupole_moments[(mol, conformer)] = self.__ext_full_quadrupole_mom(output_string)

    def __write_and_run_command(self, mol, conformer):
        """
        Writes and runs the command for GFN.

        Returns
        -------
        out_file : :class:`str`
            Name of output file with GFN-xTB results.
        """
        xyz = 'input_structure.xyz'
        out_file = 'energy.output'
        mol.write(xyz, conformer=conformer)
        # modify memory limit
        if self.mem_ulimit:
            cmd = ['ulimit -s unlimited ;']
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.gfnxtb_path)
        cmd.append(xyz)
        # set GFN Parameterization
        if self.gfn_version != '2':
            cmd.append('--gfn')
            cmd.append(self.gfn_version)
        # turn on hessian calculation if free energy requested
        cmd.append('--hess')
        # set number of cores
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # add eletronic temp term
        if self.etemp != '300':
            cmd.append('--etemp')
            cmd.append(self.etemp)
        # write solvent section of cmd
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # write charge section of cmd
        if self.charge is not None:
            cmd.append('--chrg')
            cmd.append(self.charge)
        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # uses the shell if mem_ulimit = True and waits until
        # subproces is complete. This is required to run the mem_ulimit_cmd
        # and GFN calculation in one command, which is then closed, which
        # minimizes the risk of unrestricting the memory limits.
        sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                 shell=shell)
        f.close()
        return out_file

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of molecule `mol` using GFN-xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose energy should be claculated.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        self.total_energies[(mol, conformer)] : :class:`float`
            Total energy of system.
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
            out_file = self.__write_and_run_command(mol=mol,
                                                    conformer=conformer)
            self.__get_properties(mol=mol,
                                  conformer=conformer,
                                  output_file=out_file)
        finally:
            os.chdir(init_dir)
        return self.total_energies[(mol, conformer)]
