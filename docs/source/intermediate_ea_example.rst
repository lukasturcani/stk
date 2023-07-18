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
    because the fitness function has a very low computational cost.
    However, in many applications, the fitness function will be
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
    itself, but also by the other molecules in the population. For
    example, fitness normalization can divide the fitness value of a
    molecule by the average fitness value in the population, in order
    to get the new fitness value.

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


    def get_num_rotatable_bonds(record):
        rdkit_molecule = record.get_molecule().to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return rdkit.CalcNumRotatableBonds(rdkit_molecule)


    def get_complexity(record):
        rdkit_molecule = record.get_molecule().to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_molecule)
        return BertzCT(rdkit_molecule)


    def get_num_bad_rings(record):
        rdkit_molecule = record.get_molecule().to_rdkit_mol()
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

Our ``fitness_value`` is a :class:`tuple` of the form
``(num_rotatable_bonds, complexity, num_bad_rings)``. This is a good
start, but our fitness value must be a single number. We can achieve
this by defining a :class:`.FitnessNormalizer`.

Defining a Fitness Normalizer
=============================

Fitness normalization is a process that runs after fitness calculation.
The basic idea, is that a :class:`.FitnessCalculator` takes as a
parameter a single molecule and returns its fitness value. This
fitness value will always be the same for the same molecule.
After this, we optionally perform fitness normalization with a
:class:`.FitnessNormalizer`. The :class:`.FitnessNormalizer` takes as
a parameter the entire population and yields a new population, holding
new fitness values. When a :class:`.FitnessNormalizer` assigns a
new fitness value to a molecule, its fitness value depends both
on the molecule itself, and, if desired, the fitness values of all
other molecules in the population. For example, :class:`.DivideByMean`
is a  :class:`.FitnessNormalizer`, which assigns new fitness values
according to the formula

.. code-block::

    new_fitness_value = old_fitness_value / mean_fitness_value

In order to calculate ``mean_fitness_value``, we have to be able to
consider all the fitness values in the population.

It is quite common to want to do multiple fitness normalizations in
sequence, and for this there is :class:`.NormalizerSequence`.
:class:`.NormalizerSequence` is a compound :class:`.FitnessNormalizer`.
This means it is initialized with other fitness normalizers, and
when its :meth:`~.FitnessNormalizer.normalize` method is called,
it delegates the normalization to them. For example, if we want
a :class:`.FitnessNormalizer` that first divides fitness values by
the mean fitness value in the population and then takes the inverse
of each fitness value in the population we could define

.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.DivideByMean(),
            stk.Power(-1),
        ),
    )


For the EA in this example, we want to perform a couple of
normalization steps, recall that the initial fitness values have the
form  ``(num_rotatable_bonds, complexity, num_bad_rings)``.
First, we will use :class:`.DivideByMean`, which in cases where
the fitness value is a :class:`tuple`, divides each member of the
:class:`tuple` by its own mean. This means the number of rotatable
bonds of each molecule is divided by the mean number of rotatable
bonds in the population, the complexity of each molecule is
divided by the mean complexity in the population and so on.
After using :class:`.DivideByMean`, each fitness value
is still a :class:`tuple`, but the value for each component is scaled
by the population average. This scaling is important,
because normally the different properties of a molecule have
very different orders of magnitude, which makes them very hard to
combine reasonably into a single value. However, scaling by
the population average removes differences in orders of magnitude,
and also removes the units of each quantity. This means they can
be safely combined by something like a sum. Therefore, our
initial fitness normalizer looks like this

.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.DivideByMean(),
        ),
    )

However, you might notice an issue here. We are dividing by the
mean, but the property values we are using, such as
the number of bad rings or number of rotatable bonds have values
which are allowed to be zero. This means that it's quite possible for
the population mean to be zero. If the population mean is zero
and we divide by zero - we will have a problem. We can prevent this
by adding ``1`` to every element of the :class:`tuple` before
using :class:`.DivideByMean`


.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.Add((1, 1, 1)),
            stk.DivideByMean(),
        ),
    )

:class:`.Add` is a fitness normalizer, which adds a number to
every fitness value. The number can be a :class:`tuple` of
numbers, if a our fitness value is also a :class:`tuple`.


Next, we can multiply each component of the :class:`tuple` by a
different coefficient. This will make each component contribute a
different amount to the final fitness value

.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.Add((1, 1, 1)),
            stk.DivideByMean(),
            stk.Multiply((1, 1, 1)),
        ),
    )

In our example, we do not actually have to use :class:`.Multiply`,
because all the coefficients are set to ``1``. However, if we wanted
the number of bad rings to contribute twice as much to the final
fitness value as the other properties, we would have used

.. code-block:: python

    stk.Multiply((1, 1, 2))


Now we want to combine the elements of :class:`tuple` into a single
fitness value. We can do this by taking a sum

.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.Add((1, 1, 1)),
            stk.DivideByMean(),
            stk.Multiply((1, 1, 1)),
            stk.Sum(),
        ),
    )


:class:`.Sum` is a :class:`.FitnessNormalizer`, which assumes that
each fitness value in the population is a :class:`tuple`. It then
replaces the fitness value with the sum of all elements of the
:class:`tuple`, to get the new fitness value.

Finally, we recognize that all the elements of the :class:`tuple`,
the number of rotatable bonds, complexity and the number of bad rings,
indicate a *low* fitness value. This means our final fitness value
should be the inverse of the sum, because as these values grow bigger,
our fitness value should become smaller

.. code-block:: python

    fitness_normalizer = stk.NormalizerSequence(
        fitness_normalizers=(
            stk.Add((1, 1, 1)),
            stk.DivideByMean(),
            stk.Multiply((1, 1, 1)),
            stk.Sum(),
            stk.Power(-1),
        ),
    )


That's it, our fitness normalizer will perform these steps in
order to get the final fitness values at each generation. When a
new generation is started, the fitness values of all molecules in the
population are set to the values returned by the
:class:`.FitnessCalculator`, and the fitness normalization is
started from scratch. This means that the final fitness value
can be different at each generation, even though the
:class:`.FitnessCalculator` always returns the same value for a given
molecule. This can happen, for example, because the mean value of each
member of the fitness :class:`tuple` can change at each generation,
based on different molecules being present in different generations.


Caching Fitness Values
======================

We now return to our :class:`.FitnessCalculator`, recall

.. code-block:: python

    fitness_calculator = stk.PropertyVector(
        property_functions=(
            get_num_rotatable_bonds,
            get_complexity,
            get_num_bad_rings,
        ),
    )

One of the improvements we want to make, is to prevent the
recalculation of fitness values, for molecules where a fitness value
was already found. We can do this by modifying our functions to
retrieve values from our database

.. code-block:: python

    def get_num_rotatable_bonds(database, record):
        key = stk.Inchi().get_key(record.get_molecule())
        path = "$.ea.num_rotatable_bonds"
        num_rotatable_bonds = database.get_property(key, path)
        if num_rotatable_bonds is None:
            rdkit_molecule = record.get_molecule().to_rdkit_mol()
            rdkit.SanitizeMol(rdkit_molecule)
            num_rotatable_bonds = rdkit.CalcNumRotatableBonds(rdkit_molecule)
            database.set_property(key, path, num_rotatable_bonds)
        return num_rotatable_bonds


    def get_complexity(database, record):
        key = stk.Inchi().get_key(record.get_molecule())
        path = "$.ea.complexity"
        complexity = database.get_property(key, path)
        if complexity is None:
            rdkit_molecule = record.get_molecule().to_rdkit_mol()
            rdkit.SanitizeMol(rdkit_molecule)
            complexity = BertzCT(rdkit_molecule)
            database.set_property(key, path, complexity)
        return complexity


    def get_num_bad_rings(database, record):
        key = stk.Inchi().get_key(record.get_molecule())
        path = "$.ea.num_bad_rings"
        num_bad_rings = database.get_property(key, path)
        if num_bad_rings is None:
            rdkit_molecule = record.get_molecule().to_rdkit_mol()
            rdkit.SanitizeMol(rdkit_molecule)
            num_bad_rings = sum(
                1
                for ring in rdkit.GetSymmSSSR(rdkit_molecule)
                if len(ring) < 5
            )
            database.set_property(key, path, num_bad_rings)
        return num_bad_rings


We then bind the first argument to our previously created database
object.

.. code-block:: python

    from functools import partial

    fitness_calculator = stk.PropertyVector(
        property_functions=(
            partial(get_num_rotatable_bonds, db),
            partial(get_complexity, db),
            partial(get_num_bad_rings, db),
        ),
    )


Plotting Selection
==================

One final thing we would like to do, is check which molecules
were selected for mutation, crossover and the next generation by
our EA at each generation. We can do this by creating a
:class:`.SelectionPlotter` for each :class:`.Selector` whose
selections we want to plot.

.. code-block:: python

    generation_selector = stk.Best(
        num_batches=25,
        duplicate_molecules=False,
    )
    stk.SelectionPlotter('generation_selection', generation_selector)

    mutation_selector = stk.Roulette(
        num_batches=5,
        random_seed=generator,
    )
    stk.SelectionPlotter('mutation_selection', mutation_selector)

    crossover_selector = stk.Roulette(
        num_batches=3,
        batch_size=2,
        random_seed=generator,
    )
    stk.SelectionPlotter('crossover_selection', crossover_selector)

We don't have to assign a :class:`.SelectionPlotter` to a variable,
we just have to create the instance, and it will make plots of which
molecules were selected in each :meth:`~.Selector.select` call.
This is an example graph, showing which molecules were
selected for mutation

.. image:: https://i.imgur.com/rx8qayL.png


You will get a new graph written for every :meth:`~.Selector.select`
call, where the base name of the file is determined by the string
we provided to the :class:`.SelectionPlotter`.


Final Version
=============

When we combine the code in this tutorial with the code from the
basic tutorial, we get our final version, which is best seen here_.

.. _here: https://github.com/lukasturcani/intermediate_ea/blob/master/intermediate_ea.py

Here is a plot of how the fitness value change across generations


.. image:: https://i.imgur.com/IMoMF4n.png


and here is the change in the number of rotatable bonds

.. image:: https://i.imgur.com/pWVMB14.png


We can see that our EA was pretty good lowering the number the number
of rotatable bonds, without us having to resort to any magic numbers
in our fitness function.

Finally, we can compare the initial population

.. image:: https://i.imgur.com/0px3bL0.png

to the final population

.. image:: https://i.imgur.com/wYHHUBP.png

Once again, we can see that the complexity of the molecules in the
final generation is reduced, when compared to the initial population.
