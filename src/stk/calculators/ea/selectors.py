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
from functools import wraps
from collections import Counter


logger = logging.getLogger(__name__)


class _Batch:
    """
    Represents a batch of molecules.

    """

    def __init__(self, mols, fitness_values):
        self._mols = mols
        self._fitness = sum(fitness_values[mol] for mol in mols)
        self._identity_key = frozenset(Counter(mols).items())

    def get_fitness(self):
        return self._fitness

    def get_identity_key(self):
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


def _add_yielded_reset(select):
    """
    Make every :meth:`~Selector.select` call empty :attr:`_yielded`.

    Parameters
    ----------
    select : :class:`function`
        The :meth:`~Selector.select` method to decorate.

    Returns
    -------
    :class:`function`
        The decorated :meth:`~Selector.select` method.

    """

    @wraps(select)
    def inner(self, population):
        if self._reset_yielded:
            self._yielded_mols &= set()
            self._yielded_batches &= set()
        return select(self, population)

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

        if num_batches is None:
            num_batches = float('inf')

        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        self._batch_size = batch_size
        self._num_batches = num_batches
        self._duplicate_mols = duplicate_mols
        self._duplicate_batches = duplicate_batches

        # The previously yielded molecules.
        self._yielded_mols = set()
        self._yielded_batches = set()
        # If True, _yielded is emptied before each select() call.
        self._reset_yielded = True

    def __init_subclass__(cls, **kwargs):
        cls.select = _add_yielded_reset(cls.select)
        return super().__init_subclass__(**kwargs)

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
        :class:`._Batch`
            A batch of molecules from `population`.


        """

        batches = (
            _Batch(
                mols=mols,
                fitness_values={
                    mol: fitness_values[mol] for mol in mols
                }
            )
            for mols in it.combinations(population, self._batch_size)
        )

        if not self._duplicate_batches:
            batches = filter(self._is_unyielded_batch, batches)
        if not self._duplicate_mols:
            batches = filter(self._has_unyielded_mols, batches)
        yield from batches

    def _is_unyielded_batch(self, batch):
        return batch.get_identity_key() not in self._yielded_batches

    def _has_unyielded_mols(self, batch):
        return all(mol not in self._yielded_mols for mol in batch)

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
        :class:`tuple` of :class:`.Molecule`
            A batch of selected molecules.

        """

        batches = tuple(self._get_batches(
            population=population,
            fitness_values=self._fitness_modifier(population)
        ))
        for batch in self._select(batches):
            self._yielded_batches.add(batch.get_identity_key())
            self._yielded_mols.update(batch)
            yield tuple(batch)

    def _select(self, batches):
        """
        Apply a selection algorithm to `batches`.

        Notes
        -----
        Any batch that is yielded is automatically added to
        :attr:`_yielded_mols` and :attr:`_yielded_batches`. See the
        code of :meth:`select`, which performs this operation.

        When used in a :class:`.SelectorFunnel` or
        :class:`SelectorSequence`, `batches` will not include batches
        which have been made ineligible by being yielded through a
        previous :class:`.Selector`, if :attr:`_duplicate_mols` or
        :attr:`_duplicate_batches` is set to ``False`` on the
        current :class:`.Selector`.

        Parameters
        ----------
        batches : :class:`tuple` of :class:`._Batch`
            The batches, from which some should be selected.

        Yields
        ------
        :class:`._Batch`
            A selected batch from `batches`.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and need to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class SelectorFunnel(Selector):
    """
    Applies :class:`Selector` objects in order.

    Each :class:`Selector` in :attr:`selectors` is used until
    exhaustion. The molecules selected by each :class:`Selector`
    are passed to the next one for selection. As a result, only the
    final :class:`Selector` can yield batches of size greater than 1.

    Examples
    --------
    Use :class:`Roulette` on only the 10 molecules with the highest
    fitness.

    .. code-block:: python

        import stk

        # Make a population with 20 molecules.
        pop = stk.Population(...)

        # Create a Selector which yields the top 10 molecules according
        # to roulette selection.
        best = stk.Best(num_batches=10)
        roulette = stk.Roulette()
        elitist_roulette = stk.SelectorFunnel(best, roulette)

        # Use the selector to yield molecules.
        for selected_mol in elitist_roulette.select(pop):
            # Do something with the selected molecules.
            ...

    """

    def __init__(self, *selectors):
        """
        Initialize a :class:`SelectorFunnel` instance.

        Parameters
        ----------
        *selectors : :class:`Selector`
            The :class:`Selector` objects used to select molecules. For
            all, except the last, `num_batches` must be ``1``.

        """

        self._selectors = selectors
        self._num_batches = selectors[-1]._num_batches
        self._yielded_mols = set()
        self._yielded_batches = set()
        self._reset_yielded = True
        # Make all the selectors share the same yielded set.
        for selector in self._selectors:
            # Only the funnel will reset yielded.
            selector._reset_yielded = False
            selector._yielded_mols = self._yielded_mols
            selector._yielded_batches = self._yielded_batches

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

        *head, tail = self._selectors
        for selector in head:
            population = [mol for mol, in selector.select(population)]
        yield from tail.select(population)


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

    def __init__(self, *selectors):
        """
        Initialize a :class:`SelectorSequence` instance.

        Parameters
        ----------
        *selectors : :class:`Selector`
            The :class:`Selector` objects used to select molecules.

        """

        self._selectors = selectors
        self._num_batches = sum(
            selector._num_batches for selector in self._selectors
        )
        self._yielded_mols = set()
        self._yielded_batches = set()
        self._reset_yielded = True
        # Make all the selectors share the yielded set.
        for selector in self._selectors:
            # Only the sequence will reset yielded.
            selector._reset_yielded = False
            selector._yielded_mols = self._yielded_mols
            selector._yielded_batches = self._yielded_batches

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

        for selector in self._selectors:
            yield from selector.select(population)


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
            fitness_modifier=fitness_modifier
        )

    def _select(self, batches):
        batches = sorted(batches, reverse=True)
        if not self._duplicate_batches:
            batches = filter(self._is_unyielded_batch, batches)
        if not self._duplicate_mols:
            batches = filter(self._has_unyielded_mols, batches)
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
        num_batches=None,
        duplicates=True,
        use_rank=False,
        batch_size=1,
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

        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded in more than
            one batch.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._generator = np.random.RandomState(random_seed)
        super().__init__(
            num_batches=num_batches,
            duplicates=duplicates,
            use_rank=use_rank,
            batch_size=batch_size
        )

    def select(self, population):
        """
        Yield batches of molecules using roulette selection.

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

        if not self._duplicates:
            valid_pop = [
                mol for mol in population if mol not in self._yielded
            ]
        else:
            valid_pop = list(population)
        yields = 0
        while (
            len(valid_pop) >= self._batch_size
            and yields < self._num_batches
        ):

            if self._use_rank:
                ranks = range(1, len(valid_pop)+1)
                total = sum(1/rank for rank in ranks)

                ranks = range(1, len(valid_pop)+1)
                weights = [1/(rank*total) for rank in ranks]

                valid_pop = sorted(
                    valid_pop, key=lambda m: m.fitness, reverse=True
                )

            else:
                total = sum(mol.fitness for mol in valid_pop)
                weights = [mol.fitness / total for mol in valid_pop]

            selected = tuple(self._generator.choice(
                a=valid_pop,
                size=self._batch_size,
                replace=False,
                p=weights
            ))
            yield selected
            yields += 1
            if not self._duplicates:
                self._yielded.update(selected)
                valid_pop = [
                    mol for mol in population
                    if mol not in self._yielded
                ]


class AboveAverage(Selector):
    """
    Yields above average batches of molecules.

    The fitness of a batch is the sum of all fitness values of the
    molecules in the batch. Contrary to the name, this selector will
    also yield a batch which has exactly average fitness.

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
        duplicate_batches=False,
        num_batches=None,
        duplicates=True,
        use_rank=False,
        batch_size=1
    ):
        """
        Initialize a :class:`AboveAverage` instance.

        Parameters
        ----------
        duplicate_batches : :class:`bool`, optional
            If ``True``, the same batch may be yielded more than
            once. For example, if the batch has a fitness twice above
            average it will be yielded twice, if it has a fitness
            three times above average, it will be yielded three times
            and so on.

        num_batches : :class:`int`, optional
            The number of batches to yield. If ``None`` then yielding
            will continue forever or until the generator is exhausted,
            whichever comes first.

        duplicates : :class:`bool`, optional
            If ``True``, the same molecule can be yielded in more than
            one batch.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        """

        self._duplicate_batches = duplicate_batches
        super().__init__(
            num_batches=num_batches,
            duplicates=duplicates,
            use_rank=use_rank,
            batch_size=batch_size
        )

    def select(self, population):
        """
        Yield above average fitness batches of molecules.

        Parameters
        ----------
        population : :class:`.Population`
            The population from which individuals are to be selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            The next selected batch of molecules.

        """

        mean = np.mean([
            fitness for batch, fitness in self._batch(population)
        ])

        # Sort the batches so that highest fitness batches are
        # yielded first. This is necessary because if not all batches
        # will be yielded, we want the molecules to appear in the
        # batches of the highest fitness. The same is true for when
        # self.duplicates is False. If duplicates is False then we want
        # molecules to appear in their optimal batch only.
        batches = sorted(
            self._batch(population),
            reverse=True,
            key=lambda x: x[-1]
        )

        if self._duplicates:
            batches = self._no_duplicates(batches)

        yielded = 0
        for batch, fitness in batches:
            if fitness >= mean:
                n = fitness // mean if self._duplicate_batches else 1
                for i in range(int(n)):
                    yield batch
                    yielded += 1
                    if yielded >= self._num_batches:
                        return


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
