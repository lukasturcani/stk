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


