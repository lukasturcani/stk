"""
Defines energy calculations via the ``Energy`` class.

Extending MMEA: Supporting more energy calculations.
----------------------------------------------------
It may well be the case that more ways to calculate the energy of
molcules will need to be added. For example, if support for using some
external software to perform the calculations needs to be added.

In order to add a new function which calculates the energy, just add it
as a method to the ``Energy`` class. There really aren't any
requirements beyond this, just make sure the energy value is what the
method returns.

If the energy method does not fit neatly into a single method, feel
free to split it up. Make sure that any helper methods are
private, ie that their names start with a leading underscore. Only the
main method which the user will use should be public.

To calculate the energy of some ``Molecule`` object the only thing that
is needed is calling one of the ``Energy`` methods:

    >>>  molecule.energy.rdkit('uff')
    OUT: 16.501

Here the `Energy.rdkit()` method was used as an example and `molecule`
is just some  ``Molecule`` instance. (Each ``Molecule`` object has an
``Energy`` instance in its `energy` attribute.)

Each ``Energy`` instance also has a `values` attribute which stores the
result of previous calculations. Using the molecule from the previous
example:

    >>>  molecule.energy.values
    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501}

IE this indicates that the `rdkit` method was called with the
`forcefield` argument set to 'uff' and the result was 16.501.

Calling any method of the ``Energy`` class updates the dictionary
automatically. When adding a new method to the class, no mechanism for
updating the dictionary needs to be provided. Writing the method within
the class is enough for it to update the `values` dictionary when
called.

Sometimes a function will require a parameter which does not affect the
outcome of the energy calculation. For example, to calculate the energy
using the ``macromodel`` program, the ``Energy.macromodel()`` function
can be used:

    molecule.energy.macromodel(forcefield=16,
                              macromodel_path='path/to/macromodel/dir')

This function requires the number of a forcefield (16) and the
directory where ``macromodel`` is installed on the users computer
('path/to/macromodel/dir'). However, the directory does not affect the
value of the calculated energy. When running:

    >>> molecule.energy.values

We want the output to be:

    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501,
          FunctionData('macromolecule', forcefield=16): 200}

(Assuming we are still dealing with the same `molecule` instance from
the `rdkit` example, both calculated energies should be kept in the
dictionary.)

However if we just define the ``Energy.macromodel()`` function within
the energy class, and then run:

    molecule.energy.macromodel(forcefield=16,
                           macromodel_path='path/to/macromodel/dir')

The output of

    >>> molecule.energy.values

will be

    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501,
          FunctionData('macromolecule', forcefield=16,
                        macromodel_path='path/to/macromodel/dir'): 200}

In order to make sure that the macromodel_path is excluded from the
key, decorate the function ``Energy.macromodel()`` with the
``exclude()`` decorator. For example:

    @exclude('macromodel_path')
    def macromodel(self, forcefield, macromodel_path):
        ...

Now the parameter `macromodel_path` will not form part of the key in
the `values` dictionary. If there were 2 parameters you wanted to
exlude:

    @exclude('exclude1', 'exclude2')
    def energy_func(include1, include2, exclude1, exclude2):
        ...

and so on.

Exclusion of some parameters is doubly beneficial because it means less
typing needs to be done to make the key:

    >>>  key = FunctionData('macromodel', forcefield=16)
    >>>  molecule.energy.values[key]
    OUT: 200

rather than:

    >>>  key = FunctionData('macromodel', forcefield=16,
                             macromodel_path='path/to/macromodel/dir')
    >>>  molecule.energy.values[key]
    OUT: 200

Sometimes even this isn't enough to get the key to look exactly the way
we want. For exmple:

    @exclude('force_e_calc')
    def formation(self, func, products,
                  building_blocks=None, force_e_calc=False):

The `func` parameter of this function is a FunctionData instance, which
holds data for one of the other ``Energy`` functions. This includes
all data the function requires to run, even software directories if
needed. As a result, the key when running this function may look like
this:

    FunctionData('pseudoformation', building_blocks=None,
                 func=FunctionData('macromodel', forcefield=16,
            macromodel_path='C:\\Program Files\\Schrodinger2016-3'))

Notice that the path of the macromodel installation was kept nested in
the key of the ``Energy.formation()`` function. This is undesirable. To
make this work properly, a completely custom mechanism for making the
key of the ``Energy.formation()`` is necessary. To do this, define a
function in this module. For example:

    def formation_key(fargs, fkwargs):
        ...

And set the `key` attribute of the energy function to the newly defined
function:

    Energy.formation.key = formation_key

The `fargs` and `fkwargs` arguments are the arguments and keyword
arguments with which ``Energy.formation()`` was called (including
`self`). The ``formation_key()`` function should return a FunctionData
instance which will act as the key. In our case the function was
defined so that the key is:

    FunctionData('formation', products=[], building_blocks=None,
                 func=FunctionData('macromodel', forcefield=16))

Sidenote1
---------
Make sure to use the `values` dictionary instead of running the same
calculation repeatedly.

Sidenote2
---------
The automatic updating of the dictionary is  achieved by the ``EMeta``
metaclass, ``EMethod`` descriptor and ``e_logger`` decorator. You do
not need to worry about these, but information about how they work is
provided in the docstring of ``EMeta``.

"""

import os
import rdkit.Chem.AllChem as rdkit
import subprocess as sp
import time
import psutil
import copy
from uuid import uuid4
from types import MethodType
from functools import wraps
from inspect import signature as sig


from ..convenience_tools import FunctionData
from .optimization.mopac import mopac_opt


class _EnergyError(Exception):
    """
    A class for errors in Energy methods.

    """

    ...


class EMethod:
    """
    A descriptor for methods of the ``Energy`` class.

    Attributes
    ----------
    func : function
        The method which the descriptor acts as a getter for.

    """

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, cls):
        """
        Returns a modified `self.func`.

        Attributes
        ----------
        obj : object
            The object to which the method in `self.func` should be
            bound.

        cls : object
            The class of `obj`.

        Returns
        -------
        BoundMethod
            A decorated version of the method held in `self.func`. The
            difference is that when calling the method  now, it will
            automatically update the `values` attribute of `obj`.

        self.func : function
            If the method in `self.func` is called as a class attribute
            rather than an instance attribute, return it instead of the
            descriptor.

        """
        # If the Energy method is accessed as a class attribute return
        # the unbound method.
        if obj is None:
            return self.func

        # When trying to access the Energy method as an instance
        # attribute returned a modified version of the method. The
        # method is modified so that after the method returns a value,
        # it is stored in the `values` dictionary of the Energy
        # instance.
        return e_logger(self.func, obj)


class EMeta(type):
    """
    A metaclass for ``Energy``.

    In conjuction with the EMethod descriptor and the e_logger
    decorator this function allows methods to automatically update the
    `values` attribute of their Energy instance without explicitly
    being told to do so.

    Basically this metaclass turns all methods of the Energy class into
    descriptors of the EMethod class. These descriptors return a
    decorated version of the original method defined in the class. The
    method is decorated with the ``e_logger()`` decorator. Calling this
    decorated method makes it automatically update the `values`
    dictionary of the object which used the method.

    """

    def __new__(cls, cls_name, bases, cls_dict):
        """
        Turns all the public methods in `cls` to EMethod descriptors.

        """

        # Find all the methods defined in the class `cls` and replace
        # them with an ``Emethod`` descriptor.
        for name, func in cls_dict.items():

            # Don't replace any private methods.
            if name.startswith('_') or not callable(func):
                continue

            # Replace method with descriptor.
            cls_dict[name] = EMethod(func)

        return type.__new__(cls, cls_name, bases, cls_dict)


def e_logger(func, obj):
    """
    Turns `func` into a version which updates the `values` dictionary.

    Parameters
    ----------
    func : function

    obj : Energy
        The Energy object on which the method `func` was called.

    Returns
    -------
    MethodType
        The function `func` bound to `obj` and modified so that when
        called the results update the `values` dictionary of `obj`.

    """

    @wraps(func)
    def inner(self, *args, **kwargs):

        # First get the result of the energy calculation.
        result = func(self, *args, **kwargs)

        # Next create FunctionData object to store the values of the
        # parameters used to run that calculation.
        key = func_key(func, (self,)+args, kwargs)

        # Update the `values` dictionary with the results of the
        # calculation.
        obj.values.update({key: result})
        # Return the result.
        return result

    # Make sure that the decorated function behaves like a typical
    # method and return.
    return MethodType(inner, obj)


def func_key(func, fargs=None, fkwargs=None):
    """
    Returns the key used in the `values` dictionary for `func`.

    Parameters
    ----------
    func : function
        The function whose results are to be stored in the `values`
        dictionary.

    fargs : tuple (default = None)
        The arguments passed to `func`.

    fkwargs : dict (default = None)
        The keyword arguments passed to `func`.

    Returns
    -------
    FunctionData
        The FunctionData object representing the key when `func` is
        called with the arguments `fargs` and keyword arguments
        `fkwargs`.

    """

    if fargs is None:
        fargs = []
    if fkwargs is None:
        fkwargs = []

    # Check if the function has a `key` attribute. If it does use this
    # to get its key rather than the general purpose code written here.
    if hasattr(getattr(Energy, func.__name__, False), 'key'):
        return getattr(Energy, func.__name__).key(fargs, fkwargs)

    fsig = sig(func)
    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*fargs, **fkwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key: value.default for key, value in
              dict(fsig.parameters).items() if key not in bound.keys()}
    # Combine the two sets of parameters and get rid of the `self`
    # parameter, if present.
    bound.update(default)
    if 'self' in bound.keys():
        bound.pop('self')
    # Remove any parameters that should not form key, listed in the
    # `exclude` attribute.
    if hasattr(func, 'exclude'):
        for key in func.exclude:
            bound.pop(key)

    # Return an FunctionData object representing the function and
    # chosen parameters.
    return FunctionData(func.__name__, **bound)


def exclude(*args):
    """
    A decorator to add the `exclude` attribute to methods.

    Paremeters
    ----------
    args :  tuple of strings
        Holds the names parameters which are not to be used as part of
        the key identifying the energy calculation run.

    Returns
    -------
    function
        The function which has had the `exclude` attribute added. This
        is a list  holding the names of parameters of the function
        which are not used for identifying energy calculations.

    """

    def inner(func):
        func.exclude = args
        return func

    return inner


class Energy(metaclass=EMeta):
    """
    Handles all things related to a molecules energy.

    An instance of this class will be placed in the `energy` attribute
    of each ``Molecule`` instance.

    Attributes
    ----------
    molecule : Molecule
        The energetic information held by an instance of ``Energy``
        concerns the molecule held in this attribute.

    values : dict of FunctionData instances
        The keys in the dict code for the function and parameters which
        were used to calculate the energy. The values are the energies.

    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.values = {}

    @exclude('force_e_calc')
    def formation(self, func, products,
                  building_blocks=None, force_e_calc=False):
        """
        Calculates the formation energy.

        The formation energy is calculated under the assumption that
        the molecule in `self.molecule` is composed of the molecules in
        `building_blocks` and that during formation molecules in
        `products` are formed in addition to `self.molecule`.

        Parameters
        ----------
        func : FunctionData
            A FunctionData object which describes the method of the
            ``Energy`` class used to calculate the energies of the
            various molecules. For example:

                FunctionData('rdkit', forcefield='uff')

        products : tuple of form (float, Molecule)
            This tuple holds the molecules produced in addition to
            `self.molecule`, when  a single `self.molecule` is made.
            The ``int`` represents the number made per `self.molecule`.

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in
            `self.molecule`. This argument can be omitted when the
            formation energy of a MacroMolecule instance is being
            found, as they keep this data stored elsewhere already.

        force_e_calc : bool (default = False)
            If this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.

        Returns
        -------
        float
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

        eng = self.pseudoformation(func, building_blocks, force_e_calc)
        eng -= e_products
        return eng

    def pseudoformation(self, func,
                        building_blocks=None,  force_e_calc=False):
        """
        Calculates the formation energy, sans other products.

        This is the formation energy if the energy of the other
        products of the reaction is not taken into account.

        Parameters
        ----------
        func : FunctionData
            A FunctionData object which describes the method of the
            ``Energy`` class used to calculate the energies of the
            various molecules. For example:

                FunctionData('rdkit', forcefield='uff')

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in
            `self.molecule`. This argument can be omitted when the
            formation energy of a MacroMolecule instance is being
            found, as they keep this data stored elsewhere already.

        force_e_calc : bool (default = False)
            If the this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.

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
                molf= getattr(mol.energy, func.name)
                molf(**func.params)

            e_reactants += n * mol.energy.values[fkey]

        # Get the energy of `self.molecule`. The only product whose
        # energy matters in pseudoformation.
        e_products = (self.values[fkey] if fkey in self.values.keys()
             else getattr(self, func.name)(**func.params))

        eng = e_reactants - e_products
        return eng

    def rdkit(self, forcefield):
        """
        Uses rdkit to calculate the energy of `self.molecule`.

        Parameters
        ----------
        forcefield : str
            The name of the forcefield to be used.

        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('rdkit', `forcefield`). The value is the caculated
            energy.

        Returns
        -------
        float
            The calculated energy. Note that this value is also stored
            in the dictionary `self.values`.

        """

        if forcefield == 'uff':
            self.molecule.mol.UpdatePropertyCache()
            ff = rdkit.UFFGetMoleculeForceField(self.molecule.mol)
        if forcefield == 'mmff':
            rdkit.GetSSSR(self.molecule.mol)
            self.molecule.mol.UpdatePropertyCache()
            ff = rdkit.MMFFGetMoleculeForceField(self.molecule.mol,
                  rdkit.MMFFGetMoleculeProperties(self.molecule.mol))

        eng = ff.CalcEnergy()
        return eng

    @exclude('macromodel_path')
    def macromodel(self, forcefield, macromodel_path):
        """
        Calculates the energy of `self.molecule` using macromodel.

        Note that this requires macromodel to be installed and have a
        valid license.

        Parameters
        ----------
        forcefield : int
            The id number of the forcefield to be used by macromodel.

        macromodel_path : str
            The full path of the ``Schrodinger`` suite within the
            user's machine. For example, in a default Microsoft
            installation the folder will probably be something like
            ``C:\Program Files\Schrodinger2016-2``.

        Returns
        -------
        float
            The calculated energy.

        Raises
        ------
        _EnergyError : Exception
            This exception is raised if no energy value if found in the
            MacroModel calculation's .log file. Likely due to a
            forcefield error.

        """

        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for macromodel calculations.

        # Unique file name is generated by inserting a random int into
        # the file path.
        tmp_file = "{}.mol".format(uuid4().int)
        self.molecule.write(tmp_file)

        file_root, ext = os.path.splitext(tmp_file)
        convrt_app = os.path.join(macromodel_path, 'utilities',
                                                    'structconvert')
        convrt_cmd = [convrt_app,
                     tmp_file, file_root+'.mae']
        sp.call(convrt_cmd, stdout=sp.PIPE, stderr=sp.PIPE)

        # Create an input file and run it.
        input_script = (
        "{0}.mae\n"
        "{0}-out.maegz\n"
        " MMOD       0      1      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " FFLD{1:8}      1      0      0     1.0000     0.0000     "
        "0.0000     0.0000\n"
        " BGIN       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " READ      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " ELST      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " WRIT       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " END       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n\n"
        ).format(file_root, forcefield)

        with open(file_root+'.com', 'w') as f:
            f.write(input_script)

        cmd = [os.path.join(macromodel_path, 'bmin'),
               file_root, "-WAIT", "-LOCAL"]
        sp.call(cmd)

        # Check if the license was found. If not run the function
        # again.
        with open(file_root+'.log', 'r') as f:
            log_content = f.read()
        if ('FATAL -96: Could not check out a license for mmlibs' in
            log_content):
            return self.macromodel(forcefield, macromodel_path)

        # Read the .log file and return the energy.
        with open(file_root+'.log', 'r') as f:
            for line in f:
                if "                   Total Energy =" in line:
                    eng = float(line.split()[-2].replace("=", ""))

        try:
            return eng
        except UnboundLocalError:
            raise _EnergyError('MacroModel energy calculation failed.')

    @exclude('mopac_path')
    def mopac(self, mopac_path, settings={}):
        """
        Calculates the energy of `self.molecule` using MOPAC.

        Note that this requires MOPAC to be installed and have a
        valid license.

        Parameters
        ----------
        settings: :class:`dict`, optional
        A dictionary which maps the names of the optimization parameters to
        their values. Valid values are:

            'hamiltonian' : :class:`string`
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of the NDDO
                theory, where all the atomic and diatomic parameters were
                re-optimized - update compared to PM6.
                http://openmopac.net/PM7_accuracy/PM7_accuracy.html

                Defaults to PM7.

            'method' : :class:`string`
                This calculation consists in a single point energy calculation.

            'eps' : :class:`float`
                Sets the dielectric constant for the solvent. Presence of this
                keyword will cause the COSMO (Conductor-like Screening Model)
                method to be used to approximate the effect of a solvent model
                surrounding the molecule. Solvents with low dielectric constant
                are not likely to work well with this model.
                0 means that the dielectric constant is not included in the
                calculation.
                Defaults to 80.1 can be used to model a water environment at
                room temperature.

            'charge' : :class:`list` of :class:`floats`
                When the system being studied is an ion, the charge, n, on the
                ion must be supplied as an integer. For cations n can be 1, 2,
                3, etc.; for anions -1, -2, -3, etc.
                Defaults to 0.


            'timeout' : :class:`float`
                The amount in seconds the optimization is allowed to run before
                being terminated. The default value is 2 days =
                172,800 seconds.

        mopac_path : :class:`str`
            The full path of the MOPAC suite within the
            user's machine.

        Returns
        -------
        :class:`float`
            The calculated energy.

        """
        # Define default vals for the MOPAC input
        vals = {
                'hamiltonian': 'PM7',
                'method': 'NOOPT',
                'eps': 80.1,
                'charge': 0,
                'timeout': 172800,
                }

        vals.update(settings)

        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for mopac calculations.

        # Unique file name is generated by inserting a random int into
        # the file path.
        tmp_file = "{}.mol".format(uuid4().int)
        self.molecule.write(tmp_file)

        file_root, ext = os.path.splitext(tmp_file)

        # Generate the input file
        _create_mop(file_root, self.molecule, vals)
        # Run MOPAC
        _run_mopac(file_root, mopac_path)
        return _extract_MOPAC_en(file_root)


    @exclude('mopac_path')
    def mopac_dipole(self, mopac_path, settings={}):
        """
        Calculates the dipole moment of the `self.molecule` using MOPAC.

        Note that this requires MOPAC to be installed and have a
        valid license.
        Parameters
        ----------
        settings: :class:`dict`, optional
        A dictionary which maps the names of the optimization parameters to
        their values. Valid values are:

            'hamiltonian' : :class:`string`
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of the NDDO
                theory, where all the atomic and diatomic parameters were
                re-optimized - update compared to PM6.
                http://openmopac.net/PM7_accuracy/PM7_accuracy.html

                Defaults to PM7.

            'method' : :class:`string`
                This calculation consists in a single point energy calculation.

            'eps' : :class:`float`
                Sets the dielectric constant for the solvent. Presence of this
                keyword will cause the COSMO (Conductor-like Screening Model)
                method to be used to approximate the effect of a solvent model
                surrounding the molecule. Solvents with low dielectric constant
                are not likely to work well with this model.
                0 means that the dielectric constant is not included in the
                calculation.
                Defaults to 80.1 can be used to model a water environment at
                room temperature.

            'charge' : :class:`list` of :class:`floats`
                When the system being studied is an ion, the charge, n, on the
                ion must be supplied as an integer. For cations n can be 1, 2,
                3, etc.; for anions -1, -2, -3, etc.
                Defaults to 0.


            'timeout' : :class:`float`
                The amount in seconds the optimization is allowed to run before
                being terminated. The default value is 2 days =
                172,800 seconds.

        mopac_path : :class:`str`
            The full path of the MOPAC suite within the
            user's machine.

        Returns
        -------
        :class:`float`
            The calculated dipole.

        """
        # Define default vals for the MOPAC input
        vals = {
                'hamiltonian': 'PM7',
                'method': 'NOOPT',
                'eps': 80.1,
                'charge': 0,
                'timeout': 172800,
                }

        vals.update(settings)

        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for mopac calculations.

        # Unique file name is generated by inserting a random int into
        # the file path.
        tmp_file = "{}.mol".format(uuid4().int)
        self.molecule.write(tmp_file)

        file_root, ext = os.path.splitext(tmp_file)

        # Generate the input file
        _create_mop(file_root, self.molecule, vals)
        # Run MOPAC
        _run_mopac(file_root, mopac_path)
        return _extract_MOPAC_dipole(file_root)

    @exclude('mopac_path')
    def mopac_ea(self, mopac_path, settings={}):
        """
        Calculates the electron affinity of `self.molecule` using MOPAC.
        Electron Affinity (EA): difference between energy of structure with
                                N electrons and N+1 electrons - charge = -1
                                (eV)
        Note that this requires MOPAC to be installed and have a
        valid license.

        Parameters
        ----------
        settings: :class:`dict`, optional
        A dictionary which maps the names of the optimization parameters to
        their values. Valid values are:

            'hamiltonian' : :class:`string`
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of the NDDO
                theory, where all the atomic and diatomic parameters were
                re-optimized - update compared to PM6.
                http://openmopac.net/PM7_accuracy/PM7_accuracy.html

                Defaults to PM7.

            'method' : :class:`string`
                This calculation consists in a single point energy calculation.

            'eps' : :class:`float`
                Sets the dielectric constant for the solvent. Presence of this
                keyword will cause the COSMO (Conductor-like Screening Model)
                method to be used to approximate the effect of a solvent model
                surrounding the molecule. Solvents with low dielectric constant
                are not likely to work well with this model.
                0 means that the dielectric constant is not included in the
                calculation.
                Defaults to 80.1 can be used to model a water environment at
                room temperature.

            'charge' : :class:`list` of :class:`floats`
                When the system being studied is an ion, the charge, n, on the
                ion must be supplied as an integer. For cations n can be 1, 2,
                3, etc.; for anions -1, -2, -3, etc.
                Defaults to 0.

            'timeout' : :class:`float`
                The amount in seconds the optimization is allowed to run before
                being terminated. The default value is 2 days =
                172,800 seconds.

        mopac_path : :class:`str`
            The full path of the MOPAC suite within the
            user's machine.

        mopac_path : str
            The full path of the MOPAC suite within the
            user's machine.

        Returns
        -------
        :class:`float`
            The calculated energy.

        """
        # Define default vals for the MOPAC input
        vals = {
                'hamiltonian': 'PM7',
                'method': 'NOOPT',
                'eps': 80.1,
                'charge': 0,
                'timeout': 172800,
                }
        vals.update(settings)

        # First check the energy of the neutral system
        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for mopac calculations.

        # Unique file name is generated by inserting a random int into
        # the file path.
        tmp_file = "{}.mol".format(uuid4().int)
        self.molecule.write(tmp_file)

        file_root, ext = os.path.splitext(tmp_file)

        # Generate the input file
        _create_mop(file_root, self.molecule, vals)
        # Run MOPAC
        _run_mopac(file_root, mopac_path)

        # Extract the neutral energy
        en1 = _extract_MOPAC_en(file_root)

        # Update the settings for the anion optimization
        settings2 = {
                    'method': 'OPT',
                    'gradient': 0.01,
                    'charge': -1,
                    'fileout': 'PDBOUT'
                    }

        vals.update(settings2)

        # Now generate a new molecule
        mol2 = copy.deepcopy(self.molecule)
        # Run the mopac optimisation
        mopac_opt(mol2, mopac_path, vals)
        # Extract the energy by using the self.mopac method
        vals['method'] = 'NOOPT'
        del vals['gradient']
        del vals['fileout']
        en2 = mol2.energy.mopac(mopac_path, vals)
        # Calculate the EA (eV)
        return en2 - en1


    @exclude('mopac_path')
    def mopac_ip(self, mopac_path, settings={}):
        """
        Calculates the ionization potential of `self.molecule` using MOPAC.
        Ionization Potential (IP): difference between energy of structure with
                                N electrons and N-1 electrons - charge = +1
                                (eV)
        Note that this requires MOPAC to be installed and have a
        valid license.

        Parameters
        ----------
        settings: :class:`dict`, optional
        A dictionary which maps the names of the optimization parameters to
        their values. Valid values are:

            'hamiltonian' : :class:`string`
                A series of different methods can be selected:
                PM7, PM6, AM1, CIS (CISD, CISDT), MNDO, RM1, etc..

                PM7 is the latest version of the reparametrization of the NDDO
                theory, where all the atomic and diatomic parameters were
                re-optimized - update compared to PM6.
                http://openmopac.net/PM7_accuracy/PM7_accuracy.html

                Defaults to PM7.

            'method' : :class:`string`
                This calculation consists in a single point energy calculation.

            'eps' : :class:`float`
                Sets the dielectric constant for the solvent. Presence of this
                keyword will cause the COSMO (Conductor-like Screening Model)
                method to be used to approximate the effect of a solvent model
                surrounding the molecule. Solvents with low dielectric constant
                are not likely to work well with this model.
                0 means that the dielectric constant is not included in the
                calculation.
                Defaults to 80.1 can be used to model a water environment at
                room temperature.

            'charge' : :class:`list` of :class:`floats`
                When the system being studied is an ion, the charge, n, on the
                ion must be supplied as an integer. For cations n can be 1, 2,
                3, etc.; for anions -1, -2, -3, etc.
                Defaults to 0.

            'timeout' : :class:`float`
                The amount in seconds the optimization is allowed to run before
                being terminated. The default value is 2 days =
                172,800 seconds.

        mopac_path : :class:`str`
            The full path of the MOPAC suite within the
            user's machine.

        mopac_path : str
            The full path of the MOPAC suite within the
            user's machine.

        Returns
        -------
        :class:`float`
            The calculated energy.

        """
        # Define default vals for the MOPAC input
        vals = {
                'hamiltonian': 'PM7',
                'method': 'NOOPT',
                'eps': 80.1,
                'charge': 0,
                'timeout': 172800,
                }
        vals.update(settings)

        # First check the energy of the neutral system
        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for mopac calculations.

        # Unique file name is generated by inserting a random int into
        # the file path.
        tmp_file = "{}.mol".format(uuid4().int)
        self.molecule.write(tmp_file)

        file_root, ext = os.path.splitext(tmp_file)

        # Generate the input file
        _create_mop(file_root, self.molecule, vals)
        # Run MOPAC
        _run_mopac(file_root, mopac_path)

        # Extract the neutral energy
        en1 = _extract_MOPAC_en(file_root)

        # Update the settings for the cation optimization
        settings2 = {
                    'method': 'OPT',
                    'gradient': 0.01,
                    'charge': 1,
                    'fileout': 'PDBOUT'
                    }

        vals.update(settings2)

        # Now generate a new molecule
        mol2 = copy.deepcopy(self.molecule)
        # Run the mopac optimisation
        mopac_opt(mol2, mopac_path, vals)
        # Extract the energy by using the self.mopac method
        vals['method'] = 'NOOPT'
        del vals['gradient']
        del vals['fileout']
        en2 = mol2.energy.mopac(mopac_path, vals)
        # Calculate the IP (eV)
        return en2 - en1

def formation_key(fargs, fkwargs):
    """
    Generates the key of the `formation()` method in the `values` dict.

    Parameters
    ----------
    fargs : tuple
        The arguments with which Energy.formation() was called.

    fkwargs : dict
        The keyword arguments with which Energy.formation() was called.

    Returns
    -------
    FunctionData
        The FunctionData object representing the key when `formation()`
        is called with the arguments `fargs` and keyword arguments
        `fkwargs`.

    """

    fsig = sig(Energy.formation)

    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*fargs, **fkwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key : value.default for key,value in
           dict(fsig.parameters).items() if key not in bound.keys()}

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
    Generates key of the `pseudoformation()` method in `values` dict.

    Parameters
    ----------
    fargs : tuple
        The arguments with which Energy.pseudoformation() was called.

    fkwargs : dict
        The keyword arguments with which Energy.pseudoformation() was
        called.

    Returns
    -------
    FunctionData
        The FunctionData object representing the key when
        `pseudoformation()` is called with the arguments `fargs` and
        keyword arguments `fkwargs`.

    """

    fsig = sig(Energy.pseudoformation)

    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*fargs, **fkwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key : value.default for key,value in
           dict(fsig.parameters).items() if key not in bound.keys()}

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
    To kill a MOPAC run for a specific structure it is enough to generate
    a non empty file with the molecule's name with the `.end` extension.
    """
    end_file = file_root + '.end'

    with open(end_file, 'w') as end:
        end.write('SHUT')

def _mop_line(settings):
    """
    Formats the settings dictionary with the correct keywords for MOPAC into
    a string to be added to the MOPAC input.

    Parameters
    ----------
    settings : dict
        Dictionary defined in the mopac_opt function, where all the run details
        are defined.

    Returns
    -------
    mopac_run_str : str
        String containing all the MOPAC keywords correctly formatted for the
        input file.
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
    macro_mol : MacroMolecule
        The macromolecule which is to be optimized. Its molecular
        structure file is converted to a ``.mop`` file. The original
        file is also kept.

    mopac_run_str : str
        This string specifies the MOPAC keywords to be used in the input for
        the calculation.

    Modifies
    --------
    This function creates a new ``.mop`` file from the structure file
    in `macro_mol._file`. This new file is placed in the same
    folder as the original file and has the same name with the _charge info.

    Returns
    -------
    str
        The full path of the newly created ``.mop`` file.
    """

    mop_file = file_root + '.mop'
    mol = molecule.mol

    logger.info('Creating .mop file - {}.'.format(file_root))

    # Generate the mop file containing the MOPAC run info
    with open(mop_file, 'w') as mop:
        # line for the run info
        mop.write(_mop_line(settings) + "\n")
        # line with the name of the molecule
        mop.write(file_root + "\n\n")

        # print the structural info
        for atom in mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_symbol = atom.GetSymbol()
            x, y, z = mol.GetConformer().GetAtomPosition(atom_id)
            atom_info = "{}   {}   +1  {}   +1  {}   +1 \n".format(atom_symbol,
                                                                   x, y, z)
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
