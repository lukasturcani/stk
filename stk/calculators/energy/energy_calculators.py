"""
Defines energy calculators.

.. _`adding energy calculators`:

Extending stk: Adding energy calculators.
-----------------------------------------

It may well be the case that more ways to calculate the energy of
molcules will need to be added. For example, if support for using some
external software to perform the calculations needs to be added.

In order to add a new function which calculates the energy, just add it
as a method to :class:`Energy`. There really aren't any requirements
beyond this, just make sure the energy value is what the method
returns.

If the energy method does not fit neatly into a single method, feel
free to split it up. Make sure that any helper methods are
private, i.e. that their names start with a leading underscore. Only
the main method which the user will use should be public.

To calculate the energy of some :class:`.Molecule` object the only
thing that is needed is calling one of the :class:`Energy` methods:

    >>>  molecule.energy.rdkit('uff')
    OUT: 16.501

Here, :meth:`Energy.rdkit` was used as an example and ``molecule``
is just some :class:`.Molecule` instance. Each :class:`.Molecule`
object has an :class:`Energy` instance in its :attr:`~.Molecule.energy`
attribute.

Each :class:`Energy` instance also has a :attr:`~Energy.values`
attribute, which stores the result of previous calculations. Using the
molecule from the previous example:

    >>>  molecule.energy.values
    OUT: {FunctionData('rdkit', forcefield='uff'): 16.501}

I.E. this indicates that :meth:`Energy.rdkit` was called with the
`forcefield` argument set to ``'uff'`` and the result was ``16.501``.

Calling any method of :class:`Energy` updates the dictionary
automatically. When adding a new method to the class, no mechanism for
updating the dictionary needs to be provided. Writing the method within
the class is enough for it to update the :attr:`~Energy.values`
dictionary when called.

Sometimes a function will require a parameter which does not affect the
outcome of the energy calculation. For example, to calculate the energy
using the MacroModel program, :meth:`Energy.macromodel` can be used:

.. code-block:: python

    molecule.energy.macromodel(forcefield=16,
                               macromodel_path='path/to/macromodel/dir')

This function requires the number of a forcefield (``16``) and the
directory where MacroModel is installed on the users computer
(``'path/to/macromodel/dir'``). However, the directory does not affect
the value of the calculated energy. When running:

    >>> molecule.energy.values

We want the output to be:

.. code-block:: python

    OUT: {FunctionData('rdkit', forcefield=uff): 16.501,
          FunctionData('macromodel', forcefield=16): 200}

(Assuming we are still dealing with the same ``molecule`` instance from
the :meth:`~Energy.rdkit` example, both calculated energies will be
kept in the dictionary.)

However, if we just define :meth:`Energy.macromodel` within
:class:`Energy`, and then run it:

.. code-block:: python

    molecule.energy.macromodel(forcefield=16,
                               macromodel_path='path/to/macromodel/dir')

The output of

    >>> molecule.energy.values

will be

.. code-block:: python

    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501,
          FunctionData('macromolecule',
                       forcefield=16,
                       macromodel_path='path/to/macromodel/dir'): 200}

In order to make sure that the `macromodel_path` is excluded from the
key, decorate :meth:`Energy.macromodel` with the :func:`exclude()`
decorator. For example:

.. code-block:: python

    @exclude('macromodel_path')
    def macromodel(self, forcefield, macromodel_path):
        ...

Now the parameter `macromodel_path` will not form part of the key in
the :attr:`~Energy.values` dictionary. If there were 2 parameters you
wanted to exlude:

.. code-block:: python

    @exclude('exclude1', 'exclude2')
    def energy_func(include1, include2, exclude1, exclude2):
        ...

and so on.

Sometimes even this isn't enough to get the key to look exactly the way
we want. For exmple:

.. code-block:: python

    @exclude('force_e_calc')
    def formation(self,
                  func,
                  products,
                  building_blocks=None,
                  force_e_calc=False):
        ...

The `func` parameter of this function is a :class:`.FunctionData`
instance, which holds the data of one of the other :class:`Energy`
methods. This includes all data the method requires to run, even
software directories if needed. As a result, the key when running this
function may look like this:

.. code-block:: python

    FunctionData('pseudoformation',
                 building_blocks=None,
                 func=FunctionData('macromodel',
                                   forcefield=16,
                                   macromodel_path='/opt/schrodinger2017-2'))

Notice that the path of the MacroModel installation was kept nested in
the key of :meth:`Energy.formation`. This is undesirable. To make this
work properly, a completely custom mechanism for making the
key of the :meth:`Energy.formation` is necessary. To do this, define a
function in this module. For example:

.. code-block:: python

    def formation_key(fargs, fkwargs):
        ...

And set the :attr:`key` attribute of the energy method to the newly
defined function:

.. code-block:: python

    Energy.formation.key = formation_key

The `fargs` and `fkwargs` arguments are the arguments and keyword
arguments with which :meth:`Energy.formation` was called, including
`self`. The :func:`formation_key` function should return a
:class:`.FunctionData` instance which will act as the key. In our case
the function was defined so that the key is:

.. code-block:: python

    FunctionData('formation',
                 products=[],
                 building_blocks=None,
                 func=FunctionData('macromodel', forcefield=16))

Notes.
------

#. Make sure to use the `values` dictionary instead of running the same
   calculation repeatedly.

#. The automatic updating of the dictionary is  achieved by the
   :class:`EMeta` metaclass, :func:`EMethod` descriptor and
   :func:`e_logger` decorator. You do not need to worry about these,
   but information about how they work is provided in the docstring of
   :class:`EMeta`.

"""

import os
import rdkit.Chem.AllChem as rdkit
import subprocess as sp
import psutil
import copy
from uuid import uuid4
from inspect import signature as sig
import logging

from ..utilities import FunctionData, move_generated_macromodel_files


logger = logging.getLogger(__name__)


class EnergyCalculator:
    """
    Handles all things related to a molecule's energy.

    An instance of this class will be placed in
    :attr:`.Molecule.energy` of each :class:`.Molecule` instance.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule on which all the energy calculation are performed.

    values : :class:`dict`
        A :class:`dict` of the form

        .. code-block:: python

            values = {FunctionData('method1', param='a'): 34}

        Where the :class:`.FunctionData` instance identifies an
        :class:`Energy` method which was run and the parameters which
        were used, while the key holds the result of the calculation.

    """

    def __init__(self, molecule):
        """
        Initializes a :class:`Energy` instance.

        Parameters
        ----------
        molecules : :class:`.Molecule`
            The molecule on which all the energy calculation are
            performed.

        """

        self.molecule = molecule
        self.values = {}


class FormationEnergy(EnergyCalculator):

    @exclude('force_e_calc')
    def formation(self,
                  func,
                  products,
                  building_blocks=None,
                  force_e_calc=False,
                  conformer=-1):
        """
        Calculates the formation energy.

        The formation energy is calculated under the assumption that
        the molecule in :attr:`molecule` is composed of the molecules
        in `building_blocks` and that during formation molecules in
        `products` are formed in addition to :attr:`molecule`.

        Parameters
        ----------
        func : :class:`.FunctionData`
            A :class:`.FunctionData` object which describes the method
            used to calculate the energies. For example:

            .. code-block:: python

                func  = FunctionData('rdkit', forcefield='uff')

        products : :class:`list`
            A :class:`list` of the form

            .. code-block:: python

                products = [(4, mol1), (2, mol2)]

            Where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
            of molecules produced in addition to :attr:`molecule`
            during its formation reaction. The numbers represent the
            number of each molecule made per :attr:`molecule`.

        building_blocks : :class:`list`, optional
            A :class:`list` of the form

            .. code-block:: python

                building_blocks = [(2, mol1), (3, mol2)]

            Where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
            which react to for :attr:`molecule`. The numbers represetnt
            the amount of each building block required to form 1
            :attr:`molecule`.

            This argument can be omitted when the formation energy of a
            :class:`.MacroMolecule` instance is being found, as they
            keep this data stored elsewhere already.

        force_e_calc : :class:`bool`, optional
            If this is ``True`` then all molecules in
            `building_blocks`, `products` and :attr:`molecule` will
            have their energies recalculated. Even if the energy values
            have already been found with the chosen forcefield and
            method. If ``False`` the energy is only calculated if the
            value has not already been found.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`float`
            The formation energy.

        """

        # Get the function used to calculate the energies.
        efunc = getattr(self, func.name)
        # Get the key of the function used to calculate the energies.
        fkey = func_key(efunc, None, func.params)

        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in products:
                getattr(mol.energy, func.name)(**func.params)

        # Get the total energy of the products.
        e_products = 0
        for n, mol in products:
            # If the energy has not been calculated already, calculate
            # it now.
            if fkey not in mol.energy.values.keys():
                getattr(mol.energy, func.name)(**func.params)

            e_products += n * mol.energy.values[fkey]

        eng = self.pseudoformation(func,
                                   building_blocks,
                                   force_e_calc,
                                   conformer)
        eng -= e_products
        return eng

    def pseudoformation(self,
                        func,
                        building_blocks=None,
                        force_e_calc=False,
                        conformer=-1):
        """
        Calculates the formation energy, sans other products.

        This is the formation energy if the energy of the other
        products of the reaction is not taken into account.

        Parameters
        ----------
        func : :class:`.FunctionData`
            A :class:`.FunctionData` object which describes the method
            used to calculate the energies. For example:

            .. code-block:: python

                func  = FunctionData('rdkit', forcefield='uff')

        building_blocks : :class:`list`, optional
            A :class:`list` of the form

            .. code-block:: python

                building_blocks = [(2, mol1), (3, mol2)]

            Where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
            which react to for :attr:`molecule`. The numbers represetnt
            the amount of each building block required to form 1
            :attr:`molecule`.

            This argument can be omitted when the formation energy of a
            :class:`.MacroMolecule` instance is being found, as they
            keep this data stored elsewhere already.

        force_e_calc : :class:`bool`, optional
            If this is ``True`` then all molecules in
            `building_blocks`, `products` and :attr:`molecule` will
            have their energies recalculated. Even if the energy values
            have already been found with the chosen forcefield and
            method. If ``False`` the energy is only calculated if the
            value has not already been found.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        float : :class:`float`
            The pseudoformation energy.

        """

        # Get the function used to calculate the energies.
        efunc = getattr(self, func.name)
        # Get the key of the function used to calculate the energies.
        fkey = func_key(efunc, None, func.params)

        if building_blocks is None:
            building_blocks = ((n, mol) for mol, n in
                               self.molecule.bb_counter.items())

        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in building_blocks:
                molf = getattr(mol.energy, func.name)
                molf(**func.params)

            getattr(self, func.name)(**func.params)

        # Sum the energy of building blocks under the chosen
        # forcefield.
        e_reactants = 0
        for n, mol in building_blocks:
            # If the calculation has not been done already, perform it.
            if fkey not in mol.energy.values.keys():
                molf = getattr(mol.energy, func.name)
                molf(**func.params)

            e_reactants += n * mol.energy.values[fkey]

        # Get the energy of `self.molecule`. The only product whose
        # energy matters in pseudoformation.
        e_products = (self.values[fkey] if fkey in self.values.keys()
                      else getattr(self, func.name)(**func.params))

        eng = e_reactants - e_products
        return eng


class RDKitEnergy(Energy):
    def rdkit(self, forcefield, conformer=-1):
        """
        Uses ``rdkit`` to calculate the energy.

        Parameters
        ----------
        forcefield : :class:`str`
            The name of the forcefield to be used.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`float`
            The calculated energy.

        """

        logger.debug('Starting rdkit energy calculation.')
        if forcefield == 'uff':
            self.molecule.mol.UpdatePropertyCache()
            ff = rdkit.UFFGetMoleculeForceField(self.molecule.mol,
                                                confId=conformer)
        if forcefield == 'mmff':
            rdkit.GetSSSR(self.molecule.mol)
            self.molecule.mol.UpdatePropertyCache()
            ff = rdkit.MMFFGetMoleculeForceField(
                  self.molecule.mol,
                  rdkit.MMFFGetMoleculeProperties(self.molecule.mol),
                  confId=conformer)

        eng = ff.CalcEnergy()
        return eng


def formation_key(fargs, fkwargs):
    """
    Generates the key of :meth:`Energy.formation`.

    Parameters
    ----------
    fargs : :class:`tuple`
        The arguments with which :meth:`Energy.formation` was called.

    fkwargs : :class:`dict`
        The keyword arguments with which :meth:`Energy.formation` was
        called.

    Returns
    -------
    :class:`.FunctionData`
        The :class:`.FunctionData` object representing the key of
        :meth:`Energy.formation`, when called with the arguments
        `fargs` and keyword arguments `fkwargs`.

    """

    fsig = sig(Energy.formation)

    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*fargs, **fkwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key: value.default for key, value in
               dict(fsig.parameters).items() if
               key not in bound.keys()}

    # Combine the two sets of parameters and get rid of the `self`
    # parameter, if present.
    bound.update(default)
    if 'self' in bound.keys():
        bound.pop('self')

    # Replace the energy function to be used with the key of the
    # energy function to be used.
    efuncdata = bound['func']
    efunc = getattr(Energy, efuncdata.name)
    bound['func'] = func_key(efunc, None, efuncdata.params)

    # Don't want this paramter in the key as it doenst affect the
    # result.
    bound.pop('force_e_calc')

    # Return an FunctionData object representing the function and
    # chosen parameters.
    return FunctionData('formation', **bound)


Energy.formation.key = formation_key


def pseudoformation_key(fargs, fkwargs):
    """
    Generates key of the :meth:`Energy.pseudoformation`.

    Parameters
    ----------
    fargs : :class:`tuple`
        The arguments with which :meth:`Energy.pseudoformation` was
        called.

    fkwargs : :class:`dict`
        The keyword arguments with which :meth:`Energy.pseudoformation`
        was called.

    Returns
    -------
    :class:`.FunctionData`
        The :class:`.FunctionData` object representing the key when
        :meth:`Energy.pseudoformation` is called with the arguments
        `fargs` and keyword arguments `fkwargs`.

    """

    fsig = sig(Energy.pseudoformation)

    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*fargs, **fkwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key: value.default for key, value in
               dict(fsig.parameters).items() if
               key not in bound.keys()}

    # Combine the two sets of parameters and get rid of the `self`
    # parameter, if present.
    bound.update(default)
    if 'self' in bound.keys():
        bound.pop('self')

    # Replace the energy function to be used with the key of the
    # energy function to be used.
    efuncdata = bound['func']
    efunc = getattr(Energy, efuncdata.name)
    bound['func'] = func_key(efunc, None, efuncdata.params)

    # Don't want this paramter in the key as it doenst affect the
    # result.
    bound.pop('force_e_calc')

    # Return an FunctionData object representing the function and
    # chosen parameters.
    return FunctionData('pseudoformation', **bound)


Energy.pseudoformation.key = pseudoformation_key


def _run_mopac(file_root, mopac_path, timeout=3600):

    mop_file = file_root + '.mop'

    logger.info(f'Running MOPAC - {file_root}.')

    # To run MOPAC a command is issued to the console via
    # ``subprocess.Popen``. The command is the full path of the
    # ``mopac`` program.
    file_root, ext = os.path.splitext(mop_file)
    opt_cmd = [mopac_path, file_root]
    opt_proc = psutil.Popen(opt_cmd, stdout=sp.PIPE,
                            stderr=sp.STDOUT,
                            universal_newlines=True)

    try:
        if timeout:
            proc_out, _ = opt_proc.communicate(timeout=timeout)
        else:
            proc_out, _ = opt_proc.communicate()
    except sp.TimeoutExpired:
        logger.info(('Minimization took too long and was terminated '
                     'by force - {}').format(file_root))
        _kill_mopac(file_root)


def _kill_mopac(file_root):
    """
    Stops an in-progress MOPAC run.

    To kill a MOPAC run, a file with the molecule's name and a ``.end``
    extension is written.

    Parameters
    ----------
    file_root : :class:`str`
        The molecule's name.

    Returns
    -------
    None : :class:`NoneType`

    """
    end_file = file_root + '.end'

    with open(end_file, 'w') as end:
        end.write('SHUT')


def _mop_line(settings):
    """
    Formats `settings` into a MOPAC input string.

    Parameters
    ----------
    settings : :class:`dict`
        Dictionary defined in :func:`mopac_opt`, where all the run
        details are defined.

    Returns
    -------
    :class:`str`
        String containing all the MOPAC keywords correctly formatted
        for the input file.

    """

    # Generate an empty string
    mopac_run_str = ""

    # Add Hamiltonian info
    mopac_run_str = mopac_run_str + settings['hamiltonian']
    # Forcing a single point calculation
    mopac_run_str = mopac_run_str + ' NOOPT '
    # Add EPS info
    eps_info = ' EPS={} '.format(settings['eps'])
    mopac_run_str = mopac_run_str + eps_info
    # Add Charge info
    charge_info = ' CHARGE={} '.format(settings['charge'])
    mopac_run_str = mopac_run_str + charge_info

    return mopac_run_str


def _create_mop(file_root, molecule, settings):
    """
    Creates the ``.mop`` file holding the molecule to be optimized.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule which is to be optimized. Its molecular
        structure file is converted to a ``.mop`` file. The original
        file is also kept.

    settings : :class:`dict`
        Dictionary defined by the MOPAC methods, where all the run
        details are defined.

    Returns
    -------
    :class:`str`
        The full path of the newly created ``.mop`` file.

    """

    mop_file = file_root + '.mop'
    mol = molecule.mol

    logger.info('Creating .mop file - {}.'.format(file_root))

    # Generate the mop file containing the MOPAC run info.
    with open(mop_file, 'w') as mop:
        # Line for the run info.
        mop.write(_mop_line(settings) + "\n")
        # Line with the name of the molecule.
        mop.write(file_root + "\n\n")

        # Write the structural info.
        for atom in mol.GetAtoms():
            atom_id = atom.GetIdx()
            symbol = atom.GetSymbol()
            x, y, z = mol.GetConformer().GetAtomPosition(atom_id)
            atom_info = f'{symbol}   {x}   +1  {y}   +1  {z}   +1 \n'
            mop.write(atom_info)

    return mop_file


def _extract_MOPAC_en(file_root):
    mopac_out = file_root + '.arc'

    with open(mopac_out) as outfile:
        target = "TOTAL ENERGY"
        energy_str = str([x for x in outfile.readlines() if target in x][0])
        energy_val = float(energy_str.split()[3])

    return energy_val


def _extract_MOPAC_dipole(file_root):
    mopac_out = file_root + '.arc'

    with open(mopac_out) as outfile:
        target = "DIPOLE"
        energy_str = str([x for x in outfile.readlines() if target in x][0])
        energy_val = float(energy_str.split()[2])

    return energy_val
