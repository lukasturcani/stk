"""
Selection
=========

#. :class:`.Best`
#. :class:`.Roulette`
#. :class:`.AboveAverage`
#. :class:`SelectorSequence`
#. :class:`.SelectorFunnel`


Selection is carried out by :class:`Selector` objects.
Selectors are objects with a :meth:`~Selector.select`
method, which is used to select molecules from a :class:`.Population`.
Examples of how :class:`Selector` classes can be used is given their
documentation, for example :class:`Best`, :class:`Roulette` or
:class:`AboveAverage`.

:class:`Selector` can be combined to generate more complex selection
processes, such as stochastic sampling, using
:class:`.SelectorSequence` or :class:`SelectorFunnel`. Examples of how
this can happen are given in the documentation of these two classes.


.. _`adding selectors`:

Making New Selectors
--------------------

When a new :class:`Selector` class is made it must inherit
:class:`Selector` and define a :meth:`~Selection.select` method.
:meth:`~Selection.select` must take a single argument, which is the
:class:`.Population` from which molecules are selected.

"""

import itertools as it
import numpy as np
import logging
from collections import Counter


logger = logging.getLogger(__name__)


class Batch:
    """
    Represents a batch of molecules.

    Batches can be compared, the comparison is based on their
    fitness values. Batches can also be iterated through, this
    iterates through all the molecules in the batch.

    """

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


def _is_unyielded_batch(yielded_batches):
    def inner(batch):
        return batch.get_identity_key() not in yielded_batches
    return inner


def _has_unyielded_mols(yielded_mols):
    def inner(batch):
        return all(mol not in yielded_mols for mol in batch)
    return inner


class Selector:
    """
    Selects batches of molecules from a population.

    Molecules are selected in batches, and each batch is selected
    based on its fitness. The fitness of a batch is the sum of all
    fitness values of the molecules in the batch. Batches may be of
    size 1, if single molecules are to be yielded.

    """

    def __init__(
        self,
        batch_size,
        num_batches,
        duplicate_mols,
        duplicate_batches,
        fitness_modifier,
    ):
        """
        Initialize a :class:`Selector` instance.

        Parameters
        ----------
        batch_size : :class:`int`
            The number of molecules yielded at once.

        num_batches : :class:`int`
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicate_mols : :class:`bool`
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`
            Takes the population on which :meth:`select`is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None``, each molecule is mapped to the fitness value
            found in its :attr:`~.Molecule.fitness` attribute.

        """

        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        self._batch_size = batch_size
        self._num_batches = num_batches
        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches

    @staticmethod
    def _get_fitness_values(population):
        return {mol: mol.fitness for mol in population}

    def _get_batches(self, population, fitness_values):
        """
        Get batches molecules from `population`.

        Parameters
        ----------


        Yields
        ------
        :class:`.Batch`
            A batch of molecules from `population`.


        """

        yield from (
            Batch(
                mols=mols,
                fitness_values={
                    mol: fitness_values[mol] for mol in mols
                }
            )
            for mols in it.combinations(population, self._batch_size)
        )

    def select(self, population):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
        population : :class:`.Population`
            A :class:`.Population` from which batches of molecules are
            selected.

        Yields
        ------
        :class:`Batch` of :class:`.Molecule`
            A batch of selected molecules.

        """

        batches = tuple(self._get_batches(
            population=population,
            fitness_values=self._fitness_modifier(population)
        ))

        yielded_batches = set()
        yielded_mols = set()
        selected_batches = self._select(
            batches=batches,
            yileded_mols=yielded_mols,
            yielded_batches=yielded_batches,
        )
        for batch in selected_batches:
            yielded_mols.update(batch)
            yielded_batches.add(batch.get_identity_key())
            yield batch

    def _select(self, batches, yielded_mols, yielded_batches):
        """
        Apply a selection algorithm to `batches`.

        Notes
        -----
        Any batch that is yielded is automatically added to
        `_yielded_mols` and `_yielded_batches`. See the
        code of :meth:`select`, which performs this operation.

        Parameters
        ----------
        batches : :class:`tuple` of :class:`.Batch`
            The batches, from which some should be selected.

        Yields
        ------
        :class:`.Batch`
            A selected batch from `batches`.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and need to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class RemoveBatches(Selector):
    def __init__(self, remover, selector):
        self._remover = remover
        self._selector = selector

    def select(self, population):
        removed = {
            batch.get_identity_key()
            for batch in self._remover.select(population)
        }
        batches = self._selector._get_batches(
            population=population,
            fitness_values=self._selector._fitness_modifier(population)
        )
        filtered_batches = tuple(
            batch for batch in batches
            if batch.get_identity_key() not in removed
        )

        yielded_batches = set()
        yielded_mols = set()
        selected_batches = self._selector._select(
            batches=filtered_batches,
            yileded_mols=yielded_mols,
            yielded_batches=yielded_batches,
        )
        for batch in selected_batches:
            yielded_mols.update(batch)
            yielded_batches.add(batch.get_identity_key())
            yield batch


class RemoveMolecules(Selector):
    def __init__(self, remover, selector):
        self._remover = remover
        self._selector = selector

    def select(self, population):
        removed = {
            mol
            for batch in self._remover.select(population)
            for mol in batch
        }
        population = [mol for mol in population if mol not in removed]
        yield from self._selector.select(population)


class FilterBatches(Selector):
    def __init__(self, filter, selector):
        self._filter = filter
        self._selector = selector

    def select(self, population):
        valid = {
            batch.get_identity_key()
            for batch in self._filter.select(population)
        }
        batches = self._selector._get_batches(
            population=population,
            fitness_values=self._selector._fitness_modifier(population)
        )
        filtered_batches = tuple(
            batch for batch in batches
            if batch.get_identity_key() in valid
        )

        yielded_batches = set()
        yielded_mols = set()
        selected_batches = self._selector._select(
            batches=filtered_batches,
            yileded_mols=yielded_mols,
            yielded_batches=yielded_batches,
        )
        for batch in selected_batches:
            yielded_mols.update(batch)
            yielded_batches.add(batch.get_identity_key())
            yield batch


class FilterMolecules(Selector):
    def __init__(self, filter, selector):
        self._filter = filter
        self._selector = selector

    def select(self, population):
        population = [
            mol
            for batch in self._filter.select(population)
            for mol in batch
        ]
        yield from self._selector.select(population)


class SelectorSequence(Selector):
    """
    Yields from selectors in order.

    Examples
    --------
    First use :class:`Best` to select the 10 best batches and
    then use :class:`Roulette` to select batches in proportion to their
    fitness.

    .. code-block:: python

        import stk

        # Make a population.
        pop = stk.Population(...)

        # Create a Selector which yields 10 batches of molecules and
        # then uses roulette selection indefinitely.
        best = stk.Best(batch_size=3)
        roulette = stk.Roulette(batch_size=3)
        elitist_roulette = stk.SelectorSequence(best, roulette)

        # Use the selector to yield molecules.
        for selected_mol in elitist_roulette.select(pop):
            # Do something with the selected molecules.
            ...

    """

    def __init__(self, selector1, selector2, num_batches=None):
        """
        Initialize a :class:`SelectorSequence` instance.

        Parameters
        ----------
        selector1 : :class:`.Selector`

        selector2 : :class:`.Selector`

        num_batches : :class:`int`, optional

        """

        self._selector1 = selector1
        self._selector2 = selector2
        self._num_batches = num_batches

    def select(self, population):
        """
        Yield batches of molecules.

        Parameters
        ----------
        population : :class:`.Population`
            A :class:`.Population` from which batches of molecules are
            selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            A batch of selected molecules.

        """

        yield from it.islice(
            it.chain(
                self._selector1.select(population),
                self._selector2.select(population)
            ),
            self._num_batches
        )


class Best(Selector):
    """
    Selects batches of molecules, fittest first.

    The fitness of a batch is the sum of the fitness values of the
    molecules in the batch.

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
        batch_size=1,
        num_batches=None,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`Best` instance.

        Parameters
        ----------
        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select`is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None``, each molecule is mapped to the fitness value
            found in its :attr:`~.Molecule.fitness` attribute.

        """

        super().__init__(
            batch_size=batch_size,
            num_batches=num_batches,
            dupclicate_mols=duplicate_mols,
            duplicate_batches=duplicate_batches,
            fitness_modifier=fitness_modifier,
        )

    def _select(self, batches, yielded_mols, yielded_batches):
        batches = sorted(batches, reverse=True)

        if not self._duplicate_mols:
            has_unyielded_mols = _has_unyielded_mols(yielded_mols)
            batches = filter(has_unyielded_mols, batches)

        if not self._duplicate_batches:
            is_unyielded_batch = _is_unyielded_batch(yielded_batches)
            batches = filter(is_unyielded_batch, batches)

        yield from it.islice(batches, self._num_batches)


class Worst(Selector):

    def __init__(
        self,
        batch_size=1,
        num_batches=None,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`Best` instance.

        Parameters
        ----------
        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select`is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None``, each molecule is mapped to the fitness value
            found in its :attr:`~.Molecule.fitness` attribute.

        """

        super().__init__(
            batch_size=batch_size,
            num_batches=num_batches,
            dupclicate_mols=duplicate_mols,
            duplicate_batches=duplicate_batches,
            fitness_modifier=fitness_modifier
        )

    def _select(self, batches, yielded_mols, yielded_batches):
        batches = sorted(batches)

        if not self._duplicate_mols:
            has_unyielded_mols = _has_unyielded_mols(yielded_mols)
            batches = filter(has_unyielded_mols, batches)

        if not self._duplicate_batches:
            is_unyielded_batch = _is_unyielded_batch(yielded_batches)
            batches = filter(is_unyielded_batch, batches)

        yield from it.islice(batches, self._num_batches)


class Roulette(Selector):
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
        batch_size=1,
        num_batches=None,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
        random_seed=None
    ):
        """
        Initialize a :class:`Roulette` instance.

        Parameters
        ----------
        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select`is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None``, each molecule is mapped to the fitness value
            found in its :attr:`~.Molecule.fitness` attribute.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._generator = np.random.RandomState(random_seed)
        super().__init__(
            batch_size=batch_size,
            num_batches=num_batches,
            dupclicate_mols=duplicate_mols,
            duplicate_batches=duplicate_batches,
            fitness_modifier=fitness_modifier,
        )

    def _select(self, batches, yielded_mols, yielded_batches):

        has_unyielded_mols = _has_unyielded_mols(yielded_mols)
        is_unyielded_batch = _is_unyielded_batch(yielded_batches)

        while (
            len(batches) >= self._batch_size
            and len(yielded_batches) < self._num_batches
        ):

            if not self._duplicate_mols:
                batches = filter(has_unyielded_mols, batches)
            if not self._duplicate_batches:
                batches = filter(is_unyielded_batch, batches)
            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)

            total = sum(batch.get_fitness() for batch in batches)
            weights = [
                batch.get_fitness() / total for batch in batches
            ]
            yield self._generator.choice(batches, p=weights)


class AboveAverage(Selector):
    """
    Yields above average batches of molecules.

    The fitness of a batch is the sum of all fitness values of the
    molecules in the batch.

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
        batch_size=1,
        num_batches=None,
        duplicate_mols=True,
        duplicate_batches=True,
        fitness_modifier=None,
    ):
        """
        Initialize a :class:`AboveAverage` instance.

        Parameters
        ----------
        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicate_mols : :class:`bool`, optional
            If ``True`` the same molecule can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` the same batch can be yielded more than once.

        fitness_modifier : :class:`callable`, optional
            Takes the population on which :meth:`select`is called and
            returns a :class:`dict` mapping molecules in the population
            to the fitness values the :class:`.Selector` should use.
            If ``None``, each molecule is mapped to the fitness value
            found in its :attr:`~.Molecule.fitness` attribute.

        """

        super().__init__(
            batch_size=batch_size,
            num_batches=num_batches,
            dupclicate_mols=duplicate_mols,
            duplicate_batches=duplicate_batches,
            fitness_modifier=fitness_modifier,
        )

    def _select(self, batches, yielded_mols, yielded_batches):

        if not self._duplicate_mols:
            has_unyielded_mols = _has_unyielded_mols(yielded_mols)
            batches = filter(has_unyielded_mols, batches)
        if not self._duplicate_batches:
            batches = filter(is_unyielded_batch, batches)

        is_unyielded_batch = _is_unyielded_batch(yielded_batches)

        mean = np.mean([batch.get_fitness() for batch in batches])
        batches = sorted(batches)
        while (
            batches
            and batches[-1].get_fitness() > mean
            and len(yielded_batches) < self._num_batches
        ):

            batch = batches.pop()
            fitness = batch.get_fitness()
            n = int(fitness // mean) if self._duplicate_batches else 1
            for i in range(n):
                yield batch



            if not self._duplicate_mols or not self._duplicate_batches:
                batches = tuple(batches)


class Tournament(Selector):
    """
    Yields molecules by tournament selection.

    In tournament selection, a random number of members is chosen from
    the population undergo a competition. In each competition, a random
    number of batches is chosen to compete and batches with the highest
    fitness are yielded. Competitions are repeated until the total
    number of batches yielded is equal to the number of batches.
    If batch size is greater than 1, the compared fitness is the
    average fitness of the batch. This process is repeated until the
    number of yielded batches is equal to `num_batches`.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        tournament_selection = stk.Tournament(
            num_batches=5,
            batch_size=1
        )

        # Select the molecules.
        for selected, in tournament_selection.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)

    """

    def __init__(
        self,
        num_batches=None,
        duplicates=True,
        use_rank=False,
        batch_size=1,
        duplicate_batches=False,
        random_seed=None
    ):
        """
        Initialize a :class:`TournamentSelection` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicates : :class:`bool`, optional
            If ``True``, the same molecule can be yielded in more than
            one batch.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``. In tournament sampling, this
            does not affect the selection process.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        duplicate_batches: :class:`bool` optional
            If ``True``, the same batch can be yielded from the
            selection process multiple times, as the batch can be
            selected to compete in a tournament multiple times.
            If ``False``, the batch will be removed from the batch
            population once it has been selected.

        random_seed : :class:`int`, optional
            The random seed to use.

        """
        self._duplicate_batches = duplicate_batches
        self._generator = np.random.RandomState(random_seed)
        super().__init__(
            num_batches=num_batches,
            duplicates=duplicates,
            use_rank=use_rank,
            batch_size=batch_size,
        )

    def select(self, population):
        """
        Yield molecules by tournament selection.

        Parameters
        ----------
        population : :class:`.Population`
            The population from which individuals are to be selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            The next selected batch of molecules.

        """

        # Sort batches by fitness.
        batches = sorted(
            self._batch(population),
            reverse=True,
            key=lambda x: x[-1]
        )
        # Ensure duplicate molecules are not in each batch.
        if not self._duplicates:
            # Ensure batches do not contain duplicates.
            batches = list(self._no_duplicates(batches))

        yielded = 0

        # If less than two members of the batch exist,
        # the tournament cannot take place.
        while yielded < self._num_batches and len(batches) >= 2:
            # Randomly select number of members to choose from
            # population.
            num_selections = self._generator.randint(2, len(batches)+1)
            # Get the indexes of all the batches to enter the
            # tournament.
            comparison_indexes = np.random.choice(
                len(batches),
                num_selections,
                replace=False
            )
            # Compare all batches, yielding the index of the batch
            # with the highest fitness.
            selected_index = max(
                comparison_indexes,
                key=lambda index: batches[index][-1]
            )
            # Add selected to yielded.
            self._yielded.update(batches[selected_index])
            yield batches[selected_index][0]
            # If duplicate batches are not allowed, remove the yielded
            #  batch from batches.
            if not self._duplicate_batches:
                batches.pop(selected_index)
            yielded += 1
        return


class StochasticUniversalSampling(Selector):
    """
    Yields molecules by stochastic universal sampling.

    Stochastic universal sampling yields molecules by sampling the population
    across evenly spaced intervals. Members with greater fitness values occupy
    a larger area and are therefore more likely to be sampled.
    This approach means weaker members of the population
    are given a greater chance to be chosen.

    References
    ----------
    https://en.wikipedia.org/wiki/Stochastic_universal_sampling

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

        import stk

        # Make a population holding some molecules.
        pop = stk.Population(...)

        # Make the selector.
        stochastic_sampling = stk.StochasticUniversalSampling(num_batches=5)

        # Select the molecules.
        for selected, in stochastic_samplng.select(pop):
            # Do stuff with each selected molecule, like apply a
            # mutation to it to generate a mutant.
            mutant = mutator.mutate(selected)
    """

    def __init__(
        self,
        num_batches=None,
        duplicates=False,
        use_rank=False,
        batch_size=1
    ):
        """
        Initialize a :class:`StochasticUniversalSampling` instance.

        Parameters
        ----------
        num_batches : :class:`int`, optional
            The number of batches to yield. Cannot be ``None``.
            The number of positions that will be sampled evenly from the
            population.

        duplicates : :class:`bool`, optional
            If ``True``, the same molecule can be yielded in more than
            one batch.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``. The distance of each molecule along
            is then given by ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        """
        super().__init__(
            num_batches=num_batches,
            duplicates=duplicates,
            use_rank=use_rank,
            batch_size=batch_size,
        )

    def select(self, population):
        """
        Yield molecules using stochastic universal sampling.

        Parameters
        ----------
        population : :class:`.Population`
            The population from which individuals are to be selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            The next selected batch of molecules.

        """
        # Get sorted batches of the population.
        batches = sorted(
            self._batch(population),
            reverse=True,
            key=lambda x: x[-1]
        )
        if self._duplicates:
            # Ensure batches do not contain duplicates.
            batches = list(self._no_duplicates(batches))

        yielded = 0
        if self._use_rank:
            # Set distances according to the rank.
            ranks = range(1, len(batches)+1)
            total = sum(1/rank for rank in ranks)
            distances = [1/(rank*total) for rank in ranks]
        else:
            # Set distances according to the fitness.
            total = sum(fitness for _, fitness in batches)
            distances = [fitness/total for _, fitness in batches]
        # Distance to add each time to the pointer.
        distance = 1/self._num_batches

        # First pointer position is random between 0
        # and first pointer.
        first_pointer = np.random.uniform(0, distance)
        pointer = first_pointer
        cmltv_distance = 0
        while yielded < self._num_batches:
            # Track distance progress.
            for i, dist in enumerate(distances):
                cmltv_distance += dist
                # If pointer value is less than the cumulative distance,
                # yield the molecule.
                if pointer < cmltv_distance:
                    selected = batches[i][0]
                    self._yielded.update(selected)
                    yield selected
                    yielded += 1
                    break
            pointer += distance
        return
