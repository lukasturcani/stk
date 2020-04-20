=======================
Intermediate EA Example
=======================


Introduction
============

This tutorial builds on the previous, basic, EA example. Here we will
introduce a couple more features which improve the effectiveness and
customizability of our EA. Namely:

Storing Fitness Values in a Database
    In the basic example, we did not store the calculated fitness or
    property values in our database, which is a shame, because these
    values are often expensive to calculate and generally useful.

Caching Fitness Values
    In the basic example, the fitness function recalculates fitness
    values, even if it had already calculated one for the molecule in a
    previous generation. In the basic example, this is not a big deal,
    because the fitness function has very low computational cost.
    However, in many applications the fitness function will be
    expensive, likely including a structure optimization step. In
    cases like this, it is extremely important that the fitness value
    is returned from a database, instead of being recalculated,
    whenever possible. This will greatly reduce the computational
    cost of our EA.

Normalizing Fitness Values
    This tutorial will introduce you to normalizing fitness values.
    This is an optional step in the EA, which is carried out after
    the fitness values of molecules have been calculated. Like fitness
    calculation, fitness normalization assigns fitness values to
    molecules. However, unlike fitness calculation, the fitness value
    assigned to a molecule is determined not just by the molecule
    itself, but also by the other molecules in population. For example,
    fitness normalization can divide the fitness value of a molecule
    by the average fitness value in the population, in order to
    get the new fitness value.

Plotting Selection
    This tutorial will introduce you to the :class:`.SelectionPlotter`,
    which plots which molecules were selected at each generation. This
    is often useful for validation, analysis and debugging.


You can get all the code associated with the tutorial by running::

    $ git clone https://github.com/lukasturcani/intermediate_ea


Defining a New Fitness Calculator
=================================

In the basic example, we created a :class:`.FitnessCalculator` by first
defining a Python function, which we then wrapped with
:class:`.FitnessFunction`. The function we defined calculated some
property values of a molecule, namely complexity, the number of small
rings and the number of rotatable bonds, and then combined them into
a single fitness value. One problem with this approach is that
the values were combined in a very ad-hoc way, which makes it
difficult to reason about the contribution of each property to the
final fitness value.

As an alternative to this approach, we can define three
Python functions, each of which returns a property of the molecule
whose fitness value we want

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit


    def get_num_rotatable_bonds(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return rdkit.CalcNumRotatableBonds(rdkit_molecule)


    def get_complexity(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return BertzCT(rdkit_molecule)


    def get_num_bad_rings(molecule):
        rdkit_molecule = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return sum(
            1
            for ring in rdkit.GetSymmSSSR(rdkit_molecule)
            if len(ring) < 5
        )


Now, instead of defining another Python function which combines
these property values into a single value, we can create a
:class:`.FitnessCalculator` which returns a :class:`tuple` of
property values as our fitness value. This is done by using a
:class:`.PropertyVector` instead of :class:`.FitnessFunction`

.. code-block:: python

    import stk

    fitness_calculator = stk.PropertyVector(
        property_functions=(
            get_num_rotatable_bonds,
            get_complexity,
            get_num_bad_rings,
        ),
    )


When we run

    .. code-block:: python

    fitness_value = fitness_calculator.get_fitness_value(some_molecule)

Our ``fitness_value`` with a class:`tuple` of the form
``(num_rotatable_bonds, complexity, num_bad_rings)``. This is a good
start, but our fitness value must be a single number. We can achieve
this by defining a :class:`.FitnessNormalizer`.


Defining a Fitness Normalizer
=============================

Fitness normalization is process that runs after fitness calculation.
The basic idea, is that a :class:`.FitnessCalculator` calculated
the fitness value
