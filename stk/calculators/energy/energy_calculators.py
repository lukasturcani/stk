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

import rdkit.Chem.AllChem as rdkit
import logging


logger = logging.getLogger(__name__)


class EnergyError(Exception):
    """
    Indicates a failed energy calculation.

    """

    ...


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

    def energy(self, mol, conformer=-1):
        raise NotImplementedError()


class FormationEnergy(EnergyCalculator):
    """
    Calculates the formation energy.

    The formation energy is calculated under the assumption that
    the molecule in :attr:`molecule` is composed of the molecules
    in `building_blocks` and that during formation molecules in
    `products` are formed in addition to :attr:`molecule`.

    Attributes
    ----------
    products : :class:`list`
        A :class:`list` of the form

        .. code-block:: python

            products = [(4, mol1), (2, mol2)]

        Where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
        of molecules produced in addition to :attr:`molecule`
        during its formation reaction. The numbers represent the
        number of each molecule made per :attr:`molecule`.


    """

    def __init__(self,
                 energy_calculator,
                 reactants,
                 other_products,
                 reactant_conformers=None,
                 other_product_conformers=None,
                 use_cache=True):
        """
        Initializes a :class:`FormationEnergy` instance.

        Parameters
        ----------

        """

        if reactant_conformers is None:
            reactant_conformers = [-1 for i in range(len(reactants))]

        if other_product_conformers is None:
            other_product_conformers = [
                -1 for i in range(len(other_product_conformers))
            ]

        self.energy_calculator = energy_calculator
        self.reactants = reactants
        self.other_products = other_products
        self.reactant_conformers = reactant_conformers
        self.other_product_conformers = other_product_conformers
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

        product_energy = self.energy_calculator.energy(mol, conformer)
        other_products = zip(self.other_products,
                             self.other_product_conformers)
        for product, conformer in other_products:
            product_energy += self.energy_calculator.energy(product,
                                                            conformer)

        reactants = zip(self.reactants, self.reactant_conformers)
        reactant_energy = sum(
            self.energy_calculator.energy(reactant, conformer)
            for reactant, conformer in reactants
        )

        return product_energy - reactant_energy


class MMFFEnergy(EnergyCalculator):
    """

    """

    def energy(self, mol, conformer=-1):
        """

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

        return self.force_field(self.molecule.mol, confId=conformer)


class UFFEnergy(EnergyCalculator):
    def energy(self, mol, conformer=-1):
        ...
