"""
Selection
=========

#. :class:`.Best`
#. :class:`.Worst`
#. :class:`.Roulette`
#. :class:`.AboveAverage`
#. :class:`.Tournament`
#. :class:`.StochasticUniversalSampling`
#. :class:`.RemoveBatches`
#. :class:`.RemoveMolecules`
#. :class:`.FilterBatches`
#. :class:`.FilterMolecules`
#. :class:`.If`
#. :class:`.TryCatch`
#. :class:`.Sequence`
#. :class:`.Random`
#. :class:`.RaisingCalculator`

Selection is carried out by :class:`.Selector` objects.
Selectors are objects with a :meth:`~.Selector.select`
method, which is used to select batches of molecules from an
:class:`.EAPopulation`. Examples of how :class:`.Selector` classes can
be used is given their documentation, for example :class:`.Best`,
:class:.`Roulette` or :class:`.AboveAverage`.

Selectors can be combined to generate more complex selection
processes. For example, let's say we want to implement elitism.
Elitism is when the best batches are guaranteed to be selected first,
before the selection algorithm is carried out. The
:class:`.Sequence` exists precisely for this reason. It takes
two selectors and yields batches from them, one after the other

.. code-block:: python

    import stk

    population = stk.EAPopulation(...)
    elite_roulette = stk.Sequence(
        stk.Best(5),
        stk.Roulette(20),
    )
    # Select with Best first and then with Roulette.
    for batch in elite_roulette.select(population):
        # Do stuff with batch.


What if you did not want Roulette to yield any batches
selected by :class:`.Best`? The :class:`.RemoveBatches` selector can be
used for this. It takes two selectors, one called a `remover` and one
called a `selector`. It first yields batches of molecules from a
population with the `remover`. It then passes the same population
to the `selector` but prevents it from yielding any batches
selected by the `remover`

.. code-block:: python

    roulette_without_elites = stk.RemoveBatches(
        remover=stk.Best(5),
        selector=stk.Roulette(20),
    )
    # Select batches, excluding the top 5.
    for batch in roulette_without_elites.select(population):
        # Do stuff with batch.

You can combine :class:`.RemoveBatches` and :class:`.Sequence`
to get a selector which yields the top 5 batches first and then,
using roulette, selects any of the other batches

.. code-block:: python

    elite_roulette2 = stk.Sequence(
        stk.Best(5),
        roulette_without_elites,
    )


The same thing can be written more explicitly

.. code-block:: python

    elite_roulette2 = stk.SelectorSequence(
        stk.Best(5),
        stk.RemoveBatches(
            remover=stk.Best(5),
            selector=stk.Roulette(20),
        ),
    )


You can also explore other combinations with :class:`.FilterBatches`,
:class:`.FilterMolecules` and :class:`.RemoveMolecules`. Examples
using these classes are given in their docstrings.


.. _`adding selectors`:

Making New Selectors
--------------------

When a new :class:`.Selector` class is made it must inherit
:class:`.Selector` and implement any virtual methods.

"""

import itertools as it
import numpy as np
import logging
from collections import Counter

from ..base_calculators import Calculator


logger = logging.getLogger(__name__)


class Batch:
    """
    Represents a batch of molecules.

    Batches can be compared, the comparison is based on their
    fitness values. Batches can also be iterated through, this
    iterates through all the molecules in the batch.

    Examples
    --------

    Sorting batches causes them to be sorted by fitness value.

    .. code-block:: python

        batches = (Batch(...), Batch(...), Batch(...))
        sorted_batches = sorted(batches)

    Comparison is also based on fitness value

    .. code-block:: python

        batch1 = Batch(...)
        batch2 = Batch(...)
        if batch1 > batch2:
            print('batch1 has a larger fitness value than batch2.')

    Batches can be iterated through to get the molecules in the
    batch

    .. code-block:: python

        batch = Batch(...)
        for mol in batch:
            # Do stuff with mol.

    """

    __slots__ = ['_mols', '_fitness', '_identity_key']

    def __init__(self, mols, fitness_values):
        """
        Initialize a :class:`.Batch`.

        Parameters
        ----------
        mols : :class:`tuple` of :class:`.Molecule`
            The molecules which are part of the batch.

        fitness_values : :class:`dict`
            Maps each molecule in `mols` to its fitness value.

        """

        self._mols = mols
        self._fitness = sum(fitness_values[mol] for mol in mols)
        self._identity_key = frozenset(Counter(mols).items())

    def get_size(self):
        """
        Get the number of molecules in the batch.

        Returns
        -------
        :class:`int`
            The number of molecules in the batch.

        """

        return len(self._mols)

    def get_fitness(self):
        """
        Get the fitness value of the batch.

        Returns
        -------
        :class:`float`
            The fitness value.

        """

        return self._fitness

    def get_identity_key(self):
        """
        Get the identity key of the batch.

        If two batches hold the same molecules, the same number of
        times, they will have the same identity key.

        Returns
        -------
        :class:`object`
            A hashable object which can be used to compare if two
            batches are the same.

        """

        return self._identity_key

    def __iter__(self):
        return iter(self._mols)

    def __eq__(self, other):
        return self._fitness == other._fitness

    def __gt__(self, other):
        return self._fitness > other._fitness

    def __ge__(self, other):
        return self._fitness >= other._fitness

    def __lt__(self, other):
        return self._fitness < other._fitness

    def __le__(self, other):
        return self._fitness <= other._fitness

    def __repr__(self):
        return f'Batch({", ".join(str(m) for m in self._mols)})'

    def __str__(self):
        return repr(self)


class _YieldedData:
    """
    Keeps track of batches yielded by :meth:`.Selector._select`.

    Each :meth:`.Selector._select` call should be paired with a
    new :class:`._YieldedData` instance, which should be updated
    each time a new batch is yielded.

    """

    __slots__ = ['_mols', '_batches', '_num']

    def __init__(self):
        # Has all molecules yielded by _select().
        self._mols = set()
        # Has the identity_key() of all batches yielded by select().
        self._batches = set()
        # Counts the total number of times _select() has yielded.
        self._num = 0

    def update(self, batch):
        """
        Update tracked data with a new `batch`.

        Parameters
        ----------
        batch : :class:`.Batch`
            A batch yielded by :meth:`.Selector._select`.

        Returns
        -------
        :class:`_YieldedData`
            The data tracker.

        """

        self._mols.update(batch)
        self._batches.add(batch.get_identity_key())
        self._num += 1
        return self

    def get_num(self):
        """
        Get the number of times :meth:`.Selector._select` has yielded.

        Returns
        -------
        :class:`int`
            The total number of times :meth:`.Selector._select` has
            yielded.

        """

        return self._num

    def is_yielded_batch(self, batch):
        """
        Check if `batch` has already been yielded.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` has already been yielded.

        """

        return batch.get_identity_key() in self._batches

    def is_unyielded_batch(self, batch):
        """
        Check if `batch` has not been yielded.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` has not been yielded.

        """

        return batch.get_identity_key() not in self._batches

    def has_yielded_mols(self, batch):
        """
        Check if `batch` contains any previously yielded molecules.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` contains any molecules which have
            previously been yielded.

        """

        return any(mol in self._mols for mol in batch)

    def has_no_yielded_mols(self, batch):
        """
        Check if `batch` consists only of unyielded molecules.

        Parameters
        ----------
        batch : :class:`.Batch`
            The batch to check.

        Returns
        -------
        :class:`bool`
            ``True`` if `batch` does not have any previously yielded
            molecules.

        """

        return all(mol not in self._mols for mol in batch)


class Selector(Calculator):
    """
    An abstract base class for selectors.

    Selectors select batches of molecules from a population.
    Each batch is selected based on its fitness. The fitness of a
    batch is the sum of all fitness values of the molecules in the
    batch. Batches may be of size 1.

    """

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None
    ):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            A collection of molecules from which batches are selected.

        included_batches : :class:`set`, optional
            The identity keys of batches which are allowed to be
            yielded, if ``None`` all batches can be yielded. If not
            ``None`` only batches `included_batches` will be yielded.

        excluded_batches : class:`set`, optional
            The identity keys of batches which are not allowed to be
            yielded. If ``None``, no batch is forbidden from being
            yielded.

        Yields
        ------
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        """

        # This method can be used to decorate _select in the future.
        yield from self._select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )

    def _select(self, population, included_batches, excluded_batches):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            A collection of molecules from which batches are selected.

        included_batches : :class:`set`
            The identity keys of batches which are allowed to be
            yielded, if ``None`` all batches can be yielded. If not
            ``None`` only batches `included_batches` will be yielded.

        excluded_batches : class:`set`
            The identity keys of batches which are not allowed to be
            yielded. If ``None``, no batch is forbidden from being
            yielded.

        Yields
        ------
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class _BatchingSelector(Selector):
    """
    Implements a part of the :class:`.Selector` interface.

    """

    @staticmethod
    def _return_fitness_values(population):
        return population.get_fitness_values()

    def _get_batches(
        self,
        population,
        fitness_values,
        included_batches,
        excluded_batches,
    ):
        """
        Get batches molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            The molecules which are to be batched.

        fitness_values : :class:`dict`
            Maps each molecule in `population` to the fitness value
            the selection algorithm should use.

        Yields
        ------
        :class:`.Batch`
            A batch of molecules from `population`.

        """

        for mols in it.combinations(population, self._batch_size):
            batch = Batch(
                mols=mols,
                fitness_values={
                    mol: fitness_values[mol] for mol in mols
                }
            )
            is_included = self._is_included(batch, included_batches)
            is_excluded = self._is_excluded(batch, excluded_batches)
            if is_included and not is_excluded:
                yield batch

    def _is_included(self, batch, included_batches):
        if included_batches is None:
            return True
        return batch.get_identity_key() in included_batches

    def _is_excluded(self, batch, excluded_batches):
        if excluded_batches is None:
            return False
        return batch.get_identity_key() in excluded_batches

    def _select(self, population, included_batches, excluded_batches):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            A collection of molecules from which batches are selected.

        included_batches : :class:`set`
            The identity keys of batches which are allowed to be
            yielded, if ``None`` all batches can be yielded. If not
            ``None`` only batches `included_batches` will be yielded.

        excluded_batches : class:`set`
            The identity keys of batches which are not allowed to be
            yielded. If ``None``, no batch is forbidden from being
            yielded.

        Yields
        ------
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        """

        batches = tuple(self._get_batches(
            population=population,
            fitness_values=self._fitness_modifier(population),
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        ))

        yielded = _YieldedData()
        for batch in self._select_from_batches(batches, yielded):
            yielded.update(batch)
            yield batch

        cls_name = self.__class__.__name__
        logger.debug(
            f'{cls_name} yielded {yielded.get_num()} batches.'
        )

        if (
            self._num_batches is not None
            and yielded.get_num() != self._num_batches
        ):
            logger.warning(
                f'{cls_name} was asked to yield '
                f'{self._num_batches} batches but yielded '
                f'{yielded.get_num()}.'
            )

    def _select_from_batches(self, batches, yielded):
        """
        Select batches.

        Parameters
        -----------
        batches : :class:`tuple` of :class:`.Batch`
            The batches from which some are selected.

        yielded : :class:`._YieldedData`
            Holds information on all yielded molecules and batches,
            updated automatically after every yield.

        Yields
        ------
        :class:`.Batch`
            A selected batch.

        """

        raise NotImplementedError()


class RemoveBatches(Selector):
    """
    Prevents a :class:`.Selector` from selecting some batches.

    Examples
    --------
    You want to use :class:`.Roulette` selection on all but the
    5 :class:`.Worst` batches

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.RemoveBatches(
            remover=stk.Worst(5),
            selector=stk.Roulette(20),
        )
        for batch in selector.select(population):
            # Do stuff with batch. It was selected with roulette
            # selection and is not one of the worst 5 batches.

    """

    def __init__(self, remover, selector):
        """
        Initialize a :class:`.RemoveBatches` instance.

        Parameters
        ----------
        remover : :class:`.Selector`
            Selects batches of molecules, which cannot be yielded by
            `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules, except those selected by
            `remover`.

        """

        self._remover = remover
        self._selector = selector

    def _select(self, population, included_batches, excluded_batches):
        removed_batches = {
            batch.get_identity_key()
            for batch in self._remover.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        if excluded_batches is not None:
            removed_batches |= excluded_batches

        yield from self._selector.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=removed_batches,
        )


class RemoveMolecules(Selector):
    """
    Prevents a :class:`.Selector` from selecting some molecules.

    Examples
    --------
    You want to prevent any of the molecules in the :class:`.Best`
    5 batches from being selected by :class:`.Roulette`.

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.RemoveMolecules(
            remover=stk.Best(num_batches=5, batch_size=3),
            selector=stk.Roulette(num_batches=20, batch_size=3),
        )

        for batch in selector.select(population):
            # Do stuff with batch. The batch is guaranteed not to
            # contain any molecules which are found in the best 5
            # batches of size 3.

    """

    def __init__(self, remover, selector):
        """
        Initialize a :class:`.RemoveMolecules` instance.

        Parameters
        ----------
        remover : :class:`.Selector`
            Selects batches molecules, any molecule selected cannot
            be selected by `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules, not containing any molecules
            selected by `remover`.

        """

        self._remover = remover
        self._selector = selector

    def _select(self, population, included_batches, excluded_batches):
        remover_batches = self._remover.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        removed = {mol for batch in remover_batches for mol in batch}
        valid_pop = population.__class__(*(
            mol for mol in population if mol not in removed
        ))
        valid_pop.set_fitness_values_from_dict(
            fitness_values=population.get_fitness_values(),
        )
        yield from self._selector.select(
            population=valid_pop,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )


class FilterBatches(Selector):
    """
    Allows a :class:`.Selector` to select only some batches.

    Examples
    --------
    You only want the :class:`.Best` 10 batches to participate in
    :class:`.Roulette`

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.FilterBatches(
            filter=stk.Best(10),
            selector=stk.Roulette(7),
        )
        for batch in selector.select(population):
            # Do stuff with batch. It is one of the 10 best batches and
            # was selected using roulette selection.

    """

    def __init__(self, filter, selector):
        """
        Initialize a :class:`.FilterBatches` instance.

        Parameters
        ----------
        filter : :class:`.Selector`
            Selects batches which can be yielded by `selector`.

        selector : :class:`.Selector`
            Selects batches, but only if they were also selected by
            `filter`.

        """

        self._filter = filter
        self._selector = selector

    def _select(self, population, included_batches, excluded_batches):
        allowed = {
            batch.get_identity_key()
            for batch in self._filter.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        yield from self._selector.select(
            population=population,
            included_batches=allowed,
            excluded_batches=excluded_batches,
        )


class FilterMolecules(Selector):
    """
    Allows a :class:`.Selector` to select only some molecules.

    Examples
    --------
    You want to use :class:`.Roulette` on the molecules which belong
    to the :class:`.Best` 5 batches of size 3

    .. code-block:: python

        import stk

        population = stk.Population(...)
        selector = stk.FilterMolecules(
            filter=stk.Best(num_batches=5, batch_size=3),
            selector=stk.Roulette(num_batches=20, batch_size=3),
        )
        for batch in selector.select(population):
            # Do stuff with batch. All the molecules in the batch
            # belong to the top 5 batches of size 3. The batch
            # was selected using roulette selection.

    """

    def __init__(self, filter, selector):
        """
        Initialize a :class:`.FilterMolecules` instance.

        Parameters
        ----------
        filter : :class:`.Selector`
            Selects molecules which can be yielded by `selector`.

        selector : :class:`.Selector`
            Selects batches of molecules. The batches can only
            contain molecules yielded by `filter`.

        """

        self._filter = filter
        self._selector = selector

    def _select(self, population, included_batches, excluded_batches):
        valid_batches = self._filter.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        valid_pop = population.__class__(*(
            mol for batch in valid_batches for mol in batch
        ))
        valid_pop.set_fitness_values_from_dict(
            fitness_values=population.get_fitness_values(),
        )
        yield from self._selector.select(
            population=valid_pop,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )


class Best(_BatchingSelector, Selector):
    """
    Selects batches of molecules, highest fitness value first.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        best = stk.Best()

        # Select the molecules.
        for selected, in best.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        best = stk.Best(batch_size=2)

        # Select the molecules.
        for selected in best.select(pop):
            # selected is a tuple of length 2, holding the selected
            # molecules. You can do stuff with the selected molecules
            # Like apply crossover operations on them.
            offspring = list(crosser.cross(*selected))

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`.Best` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.
            Duplicate batches can occur if the same molecule is found
            multiple times in a population.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        batches = sorted(batches, reverse=True)

        if not self._duplicate_mols:
            batches = filter(yielded.has_no_yielded_mols, batches)

        if not self._duplicate_batches:
            batches = filter(yielded.is_unyielded_batch, batches)

        yield from it.islice(batches, self._num_batches)


class Worst(_BatchingSelector, Selector):
    """
    Selects batches of molecules, lowest fitness value first.

    Examples
    --------
    Select the worst 5 batches of size 3

    .. code-block:: python

        import stk

        population = stk.Population(...)
        worst = stk.Worst(num_batches=5, batch_size=3)
        for batch in worst.select(population):
            # Do stuff with batch.


    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`.Worst` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.
            Duplicate batches can occur if the same molecule is found
            multiple times in a population.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        batches = sorted(batches)

        if not self._duplicate_mols:
            batches = filter(yielded.has_no_yielded_mols, batches)

        if not self._duplicate_batches:
            batches = filter(yielded.is_unyielded_batch, batches)

        yield from it.islice(batches, self._num_batches)


class Roulette(_BatchingSelector, Selector):
    """
    Uses roulette selection to select batches of molecules.

    In roulette selection the probability a batch is selected
    is given by its fitness. If the total fitness is the sum of all
    fitness values, the chance a batch is selected is given
    by::

        p = batch fitness / total fitness,

    where ``p`` is the probability of selection and the batch
    fitness is the sum of all fitness values of molecules in the
    batch [#]_.

    References
    ----------
    .. [#] http://tinyurl.com/csc3djm

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        roulette = stk.Roulette()

        # Select the molecules.
        for selected, in roulette.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover

    .. code-block:: python

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        roulette = stk.Roulette(batch_size=2)

        # Select the molecules.
        for selected in roulette.select(pop):
            # selected is a tuple of length 2, holding the selected
            # molecules. You can do stuff with the selected molecules
            # Like apply crossover operations on them.
            offspring = list(crosser.cross(*selected))

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None
    ):
        """
        Initialize a :class:`Roulette` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        if num_batches is None:
            num_batches = float('inf')

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._generator = np.random.RandomState(random_seed)
        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        while batches and yielded.get_num() < self._num_batches:
            total = sum(batch.get_fitness() for batch in batches)
            weights = [
                batch.get_fitness() / total for batch in batches
            ]
            yield self._generator.choice(batches, p=weights)

            if not self._duplicate_mols:
                batches = filter(yielded.has_no_yielded_mols, batches)
            if not self._duplicate_batches:
                batches = filter(yielded.is_unyielded_batch, batches)
            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)


class AboveAverage(_BatchingSelector, Selector):
    """
    Yields above average batches of molecules.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        above_avg = stk.AboveAverage()

        # Select the molecules.
        for selected, in above_avg.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        above_avg = stk.AboveAverage(batch_size=2)

        # Select the molecules.
        for selected in above_avg.select(pop):
            # selected is a tuple of length 2, holding the selected
            # molecules. You can do stuff with the selected molecules
            # Like apply crossover operations on them.
            offspring = list(crosser.cross(*selected))

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize an :class:`.AboveAverage` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        mean = np.mean([batch.get_fitness() for batch in batches])
        # Yield highest fitness batches first.
        batches = sorted(batches, reverse=True)
        # Yield only batches with a fitness larger than the mean.
        batches = it.takewhile(
            lambda batch: batch.get_fitness() > mean,
            batches
        )
        # Yield batches which are multiple times better than the mean
        # multiple times.
        batches = (
            batch
            for batch in batches
            for i in range(self._get_num_duplicates(batch, mean))
        )
        # If duplicate molecules are not allowed, allow only
        # batches with no yielded molecules.
        if not self._duplicate_mols:
            batches = filter(yielded.has_no_yielded_mols, batches)
        # If duplicate batches are not allowed, allow only
        # unyielded batches.
        if not self._duplicate_batches:
            batches = filter(yielded.is_unyielded_batch, batches)
        # Limit the number of yielded batches to _num_batches.
        yield from it.islice(batches, self._num_batches)

    def _get_num_duplicates(self, batch, mean):
        if self._duplicate_batches and self._duplicate_mols:
            return int(batch.get_fitness() // mean)
        return 1


class Tournament(_BatchingSelector, Selector):
    """
    Yields batches of molecules through tournament selection.

    In tournament selection, a random number of batches is chosen from
    the population undergo a competition. In each competition, the
    batch with the highest fitness value is yielded. This is repeated
    until `num_batches` are yielded.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        tournament = stk.Tournament(
            num_batches=5,
            batch_size=1
        )

        # Select the molecules.
        for selected, in tournament.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None,
    ):
        """
        Initialize a :class:`.Tournament` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._generator = np.random.RandomState(random_seed)
        if num_batches is None:
            num_batches = float('inf')

        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._batch_size = batch_size
        self._num_batches = num_batches
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        # The tournament can only take place if there is more than 1
        # batch.
        while (
            len(batches) > 1 and yielded.get_num() < self._num_batches
        ):
            tournament_size = self._generator.randint(
                low=2,
                high=len(batches)+1
            )
            competitors = self._generator.choice(
                a=batches,
                size=tournament_size,
                replace=False
            )
            yield max(competitors)

            if not self._duplicate_mols:
                batches = filter(yielded.has_no_yielded_mols, batches)
            if not self._duplicate_batches:
                batches = filter(yielded.is_unyielded_batch, batches)
            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)


class StochasticUniversalSampling(_BatchingSelector, Selector):
    """
    Yields batches of molecules through stochastic universal sampling.

    Stochastic universal sampling lays out batches along a line, with
    each batch taking up length proportional to its fitness. It
    then creates a set of evenly spaced pointers to different points
    on the line, each of which is occupied by a batch. Batches which
    are pointed to are yielded.

    This approach means weaker members of the population
    are given a greater chance to be chosen than in
    :class:`.Roulette` selection [#]_.

    References
    ----------
    .. [#] https://en.wikipedia.org/wiki/Stochastic_universal_sampling

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        stochastic_sampling = stk.StochasticUniversalSampling(5)

        # Select the molecules.
        for selected, in stochastic_sampling.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)
    """

    def __init__(
        self,
        num_batches=None,
        batch_size=1,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None,
    ):
        """
        Initialize a :class:`.StochasticUniversalSampling` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select` is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None`` then :meth:`.EAPopulation.get_fitness_values`
            is used.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        if fitness_modifier is None:
            fitness_modifier = self._return_fitness_values

        self._generator = np.random.RandomState(random_seed)
        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size
        self._fitness_modifier = fitness_modifier

    def _select_from_batches(self, batches, yielded):
        batches = sorted(batches, reverse=True)

        # SUS may need to run multiple rounds if duplicate_mols or
        # duplicate_batches is True. This is because in each round
        # you can generate multiple pointers to the same batch or to
        # batches sharings molecules. If this happens the lower fitness
        # batch will not be yielded. Instead a second round of SUS will
        # occur with any ineligible batches removed and a reduced
        # number of pointers, to account for batches yielded in the
        # previous rounds. This will repeat until the desired number
        # of batches has been yielded, or there are no more valid
        # batches.
        while batches and yielded.get_num() < self._num_batches:
            yield from self._select_with_stochastic_universal_sampling(
                batches=batches,
                yielded=yielded,
            )

            if yielded.get_num() < self._num_batches:
                if not self._duplicate_mols:
                    batches = (
                        filter(yielded.has_no_yielded_mols, batches)
                    )
                if not self._duplicate_batches:
                    batches = (
                        filter(yielded.is_unyielded_batch, batches)
                    )
                if (
                    not self._duplicate_mols
                    or not self._duplicate_batches
                ):
                    batches = tuple(batches)

    def _select_with_stochastic_universal_sampling(
        self,
        batches,
        yielded,
    ):

        total = sum(batch.get_fitness() for batch in batches)
        batch_positions = []
        batch_position = 0
        for batch in batches:
            batch_position += batch.get_fitness()/total
            batch_positions.append(batch_position)

        num_batches = min(
            self._num_batches - yielded.get_num(),
            len(batches)
        )
        pointer_distance = 1/num_batches
        pointers = []
        pointer = self._generator.uniform(0, pointer_distance)
        for i in range(num_batches):
            pointers.append(pointer)
            pointer += pointer_distance

        batch_index = 0
        for pointer in pointers:
            while pointer > batch_positions[batch_index]:
                batch_index += 1

            batch = batches[batch_index]

            if (
                not self._duplicate_mols
                and yielded.has_yielded_mols(batch)
            ):
                continue

            if (
                not self._duplicate_batches
                and yielded.is_yielded_batch(batch)
            ):
                continue

            yield batch
