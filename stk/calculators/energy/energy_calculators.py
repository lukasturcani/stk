"""
Defines energy calculators.

Energy calculators are objects which calculate the energy of molecules.
Each :class:`EnergyCalculator` is initialized with some settings and
calculates the energy of a molecule with
:meth:`~EnergyCalculator.energy`.

.. code-block:: python

    # Energy calculators work with any Molecule objects, such as
    # StructUnit, Polymer, Cage, etc.
    mol1 = StructUnit(...)
    mol2 = Cage(...)
    mol3 = Polymer(...)

    # Create the energy calculator.
    mmff = MMFFEnergy()

    # Calculate energies of various molecules.
    mol1_energy = mmff.energy(mol1)
    # We can optionally specify a conformer.
    mol2_energy = mmff.energy(mol2, conformer=1)
    mol3_energy = mmff.energy(mol3)

By default, calling :meth:`~EnergyCalculator.energy` twice on the
will calculate the energy a second time. However, we can use the
:attr:`~EnergyCalculator.use_cache` option to prevent recalculations
when the same molecule and conformer are given to the same energy
calculator a second time.

.. code-block:: python

    caching_mmff = MMFFEnergy(use_cache=True)
    # Calculate the energy the first time.
    energy1 = caching_mmff.energy(mol1)
    # The second time, the enegy is returned directly from memory, a
    # second calculation is not run.
    energy2 = caching_mmff.enegy(mol2)

.. _`adding energy calculators`:

Extending stk: Making new energy calculators.
---------------------------------------------

New energy calculators can be made by simply making a class which
inherits the :class:`EnergyCalculator` class. In addition to this,
the new class must define a :meth:`~EnergyCalculator.energy` method.
The method must take 2 arguments, a mandatory `mol` argument and an
optional `conformer` argument. The method will then calculate and
return the energy. There are no requirements regarding how it should
go about calculating the energy.

"""

import rdkit.Chem.AllChem as rdkit
import logging
import re
from functools import wraps
import subprocess as sp
import uuid
import os
import shutil


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


class GFNXTBEnergyInvalidSolventError(Exception):
    ...


class GFNXTBEnergy(EnergyCalculator):
    """
    Uses GFN-xTB to calculate energies.

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
    TO BE FILLED IN FROM BELOW WHEN THAT IS DONE.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = GFNXTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`GFNXTBEnergy` should usually be used after optimization with some
    other method. This is because GFN-xTB only uses xyz coordinates as input
    and so will not recognize the long bonds created during assembly.
    An optimizer which can minimize these bonds should be used before
    :class:`GFNXTBEnergy`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        uff = UFF()
        uff.optimize(polymer)
        gfnxtb = GFNXTBEnergy('/opt/gfnxtb/xtb')
        gfnxtb.energy(polymer)

    ADD EXAMPLES OF DIFFERENT ENERGIES

    """
    def __init__(self,
                 gfnxtb_path,
                 gfn_version='2',
                 output_dir=None,
                 free=False,
                 num_cores=1,
                 energy_type='total',
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 use_cache=False,
                 mem_ulimit=False,
                 strict=True):
        """
        Initializes a :class:`GFNXTBEnergy` instance.

        Parameters
        ----------
        gfnxtb_path : :class:`str`
            The path to the GFN-xTB or GFN2-xTB executable.

        gfn_version : :class:`str`
            Parameterization of GFN-xTB to use.
            See https://xtb-docs.readthedocs.io/en/latest/basics.html for a
            discussion.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        free : :class:`bool`, optional
            If ``True`` :meth:`optimize` will perform a numerical hessian
            calculation on the optimized structure to give Free energy also.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`
            The number of cores for GFN-xTB to use.

        energy_type : :class:`str`, optional
            Type of energy to extract.
            `total` : total energy of system
            `free` : total free energy of system. Reqiures free=True

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without constraints on
            the stacksize. If memory issues are encountered, this should be ``True``,
            however this may raise issues on clusters.

        etemp : :class:`int`, optional
            Electronic temperature to use (in K). Defaults to 300K.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            See https://xtb-docs.readthedocs.io/en/latest/gbsa.html for options.

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit solvent.
            Options:
                normal, tight, verytight, extreme
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        charge : :class:`str`, optional
            Formal molecular charge. `-` should be used to indicate sign.

        strict : :class:`bool`, optional
            Whether to use the `--strict` during optimization, which turns all
            internal GFN-xTB warnings into errors.

        """
        self.gfnxtb_path = gfnxtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.free = free
        self.energy_type = energy_type
        self.num_cores = str(num_cores)
        self.etemp = str(etemp)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            self.valid_solvent()
        self.solvent_grid = solvent_grid
        self.charge = charge
        self.mem_ulimit = mem_ulimit
        self.strict = strict
        super().__init__(use_cache=use_cache)

    def valid_solvent(self,):
        '''Check if solvent is valid for the given GFN version.

        See https://xtb-docs.readthedocs.io/en/latest/gbsa.html for discussion.

generalized born (GB) model with solvent accessable surface (SASA) model,
available solvents are acetone, acetonitrile, benzene (only GFN1-xTB),
CH2Cl2, CHCl3, CS2, DMF (only GFN2-xTB), DMSO, ether, H2O, methanol,
n-hexane (only GFN2-xTB), THF and toluene. The solvent input is not case-sensitive.
 The Gsolv reference state can be chosen as reference or bar1M (default).

        '''
        if self.gfn_version == '0':
            raise GFNXTBEnergyInvalidSolventError(
                f'No solvent valid for version: {self.gfn_version}'
            )
        elif self.gfn_version == '1':
            valid_solvents = ['acetone', 'acetonitrile', 'benzene',
                              'CH2Cl2'.lower(), 'CHCl3'.lower(), 'CS2'.lower(),
                              'DMF'.lower(), 'DMSO'.lower(), 'ether', 'H2O'.lower(),
                              'methanol', 'THF'.lower(), 'toluene']
            if self.solvent in valid_solvents:
                return True
            else:
                raise GFNXTBEnergyInvalidSolventError(
                    f'{self.solvent} is an invalid solvent for version {self.gfn_version}!'
                )
        elif self.gfn_version == '1':
            valid_solvents = ['acetone', 'acetonitrile', 'CH2Cl2'.lower(),
                              'CHCl3'.lower(), 'CS2'.lower(), 'DMF'.lower(),
                              'DMSO'.lower(), 'ether', 'H2O'.lower(), 'methanol',
                              'n-hexane'.lower(), 'THF'.lower(), 'toluene']
            if self.solvent in valid_solvents:
                return True
            else:
                raise GFNXTBEnergyInvalidSolventError(
                    f'{self.solvent} is an invalid solvent for version {self.gfn_version}!'
                )

    def get_properties(self):
        """Extract desired properties from GFN-xTB single point energy calculation.

        """
        self.properties = {}

        # energy
        self.properties['totalenergy'] = 0.0

        return self.properties

    def extract_energy(self, output_file):
        """
        Extract desired energy from GFN2-xTB output file.

        Format depends on version. Works with version 190418.


        Obtained results (in a.u.):
            - free energies (FE)
            - absolute energy (TE)
            - SCC energy (SCE)
            - G (GT)
            - H (HT)

        """
        print('CHECK OTHER ENERGY TYPES')
        if self.free is False:
            if self.energy_type in ['total']:
                raise(f'{self.energy_type} requires hessian calculation - free=True')

        for line in open(output_file, 'r'):
            # regex:
            nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
            # free energy in a.u.
            if '| TOTAL FREE ENERGY' in line and self.energy_type == 'total':
                FE_au = nums.search(line.rstrip()).group(0)
                energy_au = FE_au
                break
            if '| TOTAL ENERGY' in line and self.energy_type == 'free':
                TE_au = nums.search(line.rstrip()).group(0)
                energy_au = TE_au
                break
            if '| TOTAL ENTHALPY' in line and self.energy_type == 'enthalpy':
                HT_au = nums.search(line.rstrip()).group(0)
                energy_au = HT_au
                break
        return energy_au

    def energy(self, mol, conformer=-1):
        """
        Calculates the energy of molecule `mol` using GFN-xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        self.properties : :class:`dict`
            Dictionary containing desired properties.

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
            xyz = 'input_structure.xyz'
            out_file = 'output_info.output'
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
            if self.free is True:
                cmd.append('--hess')
                cmd.append(self.opt_level)
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
            # add strict term
            if self.strict is True:
                cmd.append('--strict')
            cmd = ' '.join(cmd)
            print(cmd)
            f = open(out_file, 'w')
            # uses the shell if mem_ulimit = True and waits until
            # subproces is complete. This is required to run the mem_ulimit_cmd
            # and GFN calculation in one command, which is then closed, which
            # minimizes the risk of unrestricting the memory limits.
            sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                     shell=shell)
            f.close()
            self.get_properties()
        finally:
            os.chdir(init_dir)
        return self.properties
