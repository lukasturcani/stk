"""
Defines selection functions via :class:`Selection`.

.. _`adding selection functions`:

Extending stk: Adding selection functions.
------------------------------------------

If a new selection operation is to be added to ``stk`` it should be
added as a method in :class:`Selection`. The only requirements are
that the first argument is `population` (excluding any `self` or `cls`
arguments).

The naming requirement exists to help users identify which arguments
are handled automatically by the GA and which they need to define in
the input file. The convention is that if the selection function takes
an argument called  `population` it does not have to be specified in
the input file.

All selection functions should be defined as generators, which yield a
member of the population. In the case of crossover selection functions,
they should yield a tuple of members. Generational selection functions
should yield each molecule at most once.

If the selection function does not fit neatly into a single function
make sure that any helper functions are private, i.e. that their names
start with a leading underscore.

"""

import itertools as it
import numpy as np
import logging


logger = logging.getLogger(__name__)


class SelectorError(Exception):
    ...


class Selector:
    """
    Selects molecules from a population.

    Molecules are selected in batches, and each batch is selected
    based on its fitness. Batches may be of size 1, if single molecules
    are to be yielded.

    Attributes
    ----------
    duplicates : :class:`bool`
        If ``True`` the same member can be yielded in more than
        one batch.

    use_rank : :class:`bool`
        When ``True`` the fitness value of an individual is
        calculated as ``f = 1/rank``.

    batch_size : :class:`int`
        The number of molecules yielded at once.

    elitism : :class:`int`
        If not ``None`` then the number determines how many of the
        most fit molecules are exclusively used for selection. For
        example, if :attr:`elitism` is ``5`` then only the top 5
        fittest molecules can be selected.

    truncation : :class:`int`
        If not ``None`` then the number determines how many of the
        least fit molecules are removed from the population before
        selection takes place.

    """

    def __init__(self,
                 duplicates,
                 use_rank,
                 batch_size,
                 elitism,
                 truncation):
        """
        Initializes a :class:`Selector` instance.

        Parameters
        ----------
        duplicates : :class:`bool`
            If ``True`` the same member can be yielded in more than
            one batch.

        use_rank : :class:`bool`
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`
            The number of molecules yielded at once.

        elitism : :class:`int`
            If not ``None`` then the number determines how many of the
            most fit molecules are exclusively used for selection. For
            example, if :attr:`elitism` is ``5`` then only the top 5
            fittest molecules can be selected.

        truncation : :class:`int`
            If not ``None`` then the number determines how many of the
            least fit molecules are removed from the population before
            selection takes place.

        Raises
        ------
        :class:`SelectorError`
            If both `elitism` and `truncation` are not ``None``.

        """

        if elitism is not None and truncation is not None:
            raise SelectorError()

        self.batch_size = batch_size
        self.elitism = elitism
        self.truncation = truncation

    def select(self, population):
        """
        Selects molecules from `population`.

        Parameters
        ----------
        population : :class:`.Population`
            A :class:`.Population` from which molecules are selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            A batch of selected molecules.

        """

        return NotImplementedError()


class SelectorPipeline(Selector):
    """
    Applies a set of :class:`Selector`s in order.

    This selector does not inherit any attributes from its base
    class.

    Attributes
    ----------
    selectors : :class:`tuple` of :class:`Selector`
        The selectors which get used to select molecules.

    """

    def __init__(self, *selectors):
        """
        Initializes a :class:`SelectorPipeline` instance.

        Parameters
        ----------
        selectors : :class:`tuple` of :class:`Selector`
            The selectors which get used to select molecules.

        """

        self.selectors = selectors

    def select(self, population):
        """
        Selects molecules from `population`.

        Parameters
        ----------
        population : :class:`.Population`
            A :class:`.Population` from which molecules are selected.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            A batch of selected molecules.

        """

        for selector in self.selectors:
            logger.info(
                f'Using {selector.__class__.__name__} for selection.'
            )
            yield from selector.select(population)


class Fittest(Selector):
    """
    Selects molecules, fittest first.

    Examples
    --------
    Yielding molecules one at a time. For example, if molecules need
    to be selected for mutation or the next generation.

    .. code-block:: python

    Yielding multiple molecules at once. For example, if molecules need
    to be selected for crossover.

    .. code-block:: python


    """

    def __init__(self,
                 duplicates=True,
                 batch_size=1,
                 elitism=None,
                 truncation=None):
        """
        Initializes a :class:`Fittest` instance.

        Parameters
        ----------
        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded in more than
            one batch.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        elitism : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            most fit molecules are exclusively used for selection. For
            example, if :attr:`elitism` is ``5`` then only the top 5
            fittest molecules can be selected.

        truncation : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            least fit molecules are removed from the population before
            selection takes place.

        """

        super(duplicates=duplicates,
              use_rank=False,
              batch_size=batch_size,
              elitism=elitism,
              truncation=truncation)

    def select(self, population):
        """
        Yields members of the population, fittest first.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population who's members need to be selected.

        Yields
        ------
        :class:`.MacroMolecule`
            The next fittest :class:`.MacroMolecule` instance held by
            `population`.

        """

        if self.elitism is not None or self.truncation is not None:
            population = sorted(population, reverse=True)

            if self.elitism is not None:
                population = population[:self.elitism]

            elif self.truncation is not None:
                population = population[:-self.truncation]

        # For each batch, sum the fitness values of all molecules in
        # in the batch.
        batches = (
            (batch, sum(m.fitness for m in batch))
            for batch in it.combinations(population, self.batch_size)
        )

        # Sort by total fitness value of each batch.
        batches = sorted(batches, reverse=True, key=lambda x: x[1])

        batches = (batch for batch, fitness in batches)

        if not self.duplicates:
            batches = self._no_duplicates(batches)

        yield from batches

    @staticmethod
    def _no_duplicates(batches):
        """
        Makes sure that no molecules is yielded in more than one batch.

        Parameters
        ----------
        batches : :class:`iterable`
            An :class:`iterable` yielding :class:`tuple` of
            :class:`.Molecule`.

        Yields
        ------
        :class:`tuple` of :class:`.Molecule`
            Batches of molecules, none of which share molecules.

        """

        seen = set()
        for batch in batches:
            if all(mol not in seen for mol in batch):
                seen.update(batch)
                yield batch


class Roulette(Selector):
    """
    Uses roulette selection to select molecules.

    Attributes
    ----------
    duplicate_batches : :class:`bool`
        If ``True`` then the same batch can be yielded more than once.


    """

    def __init__(self,
                 duplicates=True,
                 duplicate_batches=True,
                 use_rank=False,
                 batch_size=1,
                 elitism=None,
                 truncation=None):
        """
        Initializes a :class:`Roulette` instance.

        Parameters
        ----------
        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` then the same batch can be yielded more than
            once.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        elitism : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            most fit molecules are exclusively used for selection. For
            example, if :attr:`elitism` is ``5`` then only the top 5
            fittest molecules can be selected.

        truncation : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            least fit molecules are removed from the population before
            selection takes place.

        """

        self.duplicate_batches = duplicate_batches
        super(duplicates=duplicates,
              use_rank=use_rank,
              batch_size=batch_size,
              elitism=elitism,
              truncation=truncation)

    def select(self, population):
        """
        Yields individuals using roulette selection.

        In roulette selection the probability an individual is selected
        is given by its fitness. If the total fitness is the sum of all
        fitness values, the chance an individuals is selected is given
        by::

            p = individual fitness / total fitness,

        where ``p`` is the probability of selection [#]_.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        References
        ----------
        .. [#] http://tinyurl.com/csc3djm

        """

        if self.elitism is not None or self.truncation is not None:
            population = sorted(population, reverse=True)

            if self.elitism is not None:
                population = population[:self.elitism]

            elif self.truncation is not None:
                population = population[:-self.truncation]

        yielded_mols = set()
        yielded_batches = set()

        valid_pop = list(population)
        while valid_pop:
            total = sum(ind.fitness for ind in valid_pop)
            weights = [ind.fitness / total for ind in valid_pop]

            selected = np.random.choice(valid_pop, p=weights)
            yielded.add(selected) if not duplicates else None
            yield selected
            valid_pop = [ind for ind in population if ind not in yielded]


class DeterministicSampling(Selector):
    """
    Uses deterministic sampling to select molecules.

    Attributes
    ----------
    duplicate_batches : :class:`bool`
        If ``True`` then the same batch can be yielded more than once.

    Examples
    --------

    """

    def __init__(self,
                 duplicates=True,
                 duplicate_batches=True,
                 use_rank=False,
                 batch_size=1,
                 elitism=None,
                 truncation=None):
        """
        Initializes a :class:`DeterministicSampling` instance.

        Parameters
        ----------
        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` then the same batch can be yielded more than
            once.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        elitism : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            most fit molecules are exclusively used for selection. For
            example, if :attr:`elitism` is ``5`` then only the top 5
            fittest molecules can be selected.

        truncation : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            least fit molecules are removed from the population before
            selection takes place.

        """

        self.duplicate_batches = duplicate_batches
        super(duplicates=duplicates,
              use_rank=use_rank,
              batch_size=batch_size,
              elitism=elitism,
              truncation=truncation)

    def select(self, population, truncation=None):
        """
        Yields individuals using deterministic sampling.

        This algorithm can only be used for selection of mutants. It
        will return duplicates. The results of :meth:`fittest` are
        equivalent to algorithm without duplicates.

        In determnistic sampling the mean fitness value of the
        population is calculated, ``<f>``. For each individual a
        normalized fitness value is then calculated via::

            fn = f / <f>

        where ``fn`` is the normalized fitness value and ``f`` is the
        original fitness value. An individual will be selected ``n``
        times, where ``n`` is the integer value of ``fn``. After all
        individuals where ``n > 0`` are yielded, all inidividuals are
        yielded again, largest decimal part of ``fn`` first.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        """

        if self.elitism is not None or self.truncation is not None:
            population = sorted(population, reverse=True)

            if self.elitism is not None:
                population = population[:self.elitism]

            elif self.truncation is not None:
                population = population[:-self.truncation]

        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]
        mean = np.mean([mem.fitness for mem in pop])
        decimals = []
        for mem in pop:
            q, r = divmod(mem.fitness, mean)
            decimals.append((r, mem))
            for i in range(int(q)):
                yield mem

        for r, mem in sorted(decimals, reverse=True):
            yield mem


class StochasticSampling(Selector):
    """
    Uses stochastic sampling to select molecules.

    Attributes
    ----------
    duplicate_batches : :class:`bool`
        If ``True`` then the same batch can be yielded more than once.

    """

    def __init__(self,
                 duplicates=True,
                 duplicate_batches=True,
                 use_rank=False,
                 batch_size=1,
                 elitism=None,
                 truncation=None):
        """
        Initializes a :class:`StochasticSampling` instance.

        Parameters
        ----------
        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded in more than
            one batch.

        duplicate_batches : :class:`bool`, optional
            If ``True`` then the same batch can be yielded more than
            once.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as ``f = 1/rank``.

        batch_size : :class:`int`, optional
            The number of molecules yielded at once.

        elitism : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            most fit molecules are exclusively used for selection. For
            example, if :attr:`elitism` is ``5`` then only the top 5
            fittest molecules can be selected.

        truncation : :class:`int`, optional
            If not ``None`` then the number determines how many of the
            least fit molecules are removed from the population before
            selection takes place.

        """

        self.duplicate_batches = duplicate_batches
        super(duplicates=duplicates,
              use_rank=use_rank,
              batch_size=batch_size,
              elitism=elitism,
              truncation=truncation)

    def select(self, population):
        """
        Yields individuals via stochastic sampling.

        Each fitness value is used to calculate the normalized fitness
        of an individual::

            fn = f / <f>

        where ``fn`` is the normalized fitness value, ``f`` is the
        original fitness value and ``<f>`` is the mean fitness of the
        population. If ``fn`` is greater than ``1`` then the individual
        is guaranteed to be selected. If duplicates are allowed,
        individuals are guaranteed to be selected ``n`` times, where
        ``n`` is the integer part of ``fn``.

        Once all the guarnteed individuals have been yielded, the
        remaining individuals are yielded via the roulette method. The
        weights in the roulette method are based on the decimal part of
        ``fn``.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        elites : :class:`int`, optional
            The number of the fittest members which are guaranteed to
            be yielded first.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded more than
            once.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as, ``f = 1/rank``.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        """

        if self.elitism is not None or self.truncation is not None:
            population = sorted(population, reverse=True)

            if self.elitism is not None:
                population = population[:self.elitism]

            elif self.truncation is not None:
                population = population[:-self.truncation]

        yielded = set()
        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        for x in range(elites):
            yielded.add(pop[x])
            yield pop[x]

        for r, mem in enumerate(pop, 1):
            mem._ssfitness = 1 / r if use_rank else mem.fitness

        yield from self._stochastic_sampling_guaranteed(pop,
                                                        yielded,
                                                        duplicates)
        yield from self._stochastic_sampling_roulette(pop,
                                                      yielded,
                                                      duplicates)

    def _stochastic_sampling_guaranteed(self, pop,
                                        yielded, duplicates):
        """
        Yielded the members guaranteed by stochastic sampling.

        Parameters
        ----------
        pop : :class:`list` of :class:`.MacroMolecule`
            The molecules which are being selected.

        yielded : :class:`set` of :class:`.MacroMolecule`
            Holds all previously yielded molecules.

        duplicates : :class:`bool`
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        :class:`.MacroMolecule`
            The next member whose normalized fitness integer component
            is greater than 0.

        """

        mean = np.mean([mem._ssfitness for mem in pop])
        for mem in pop:
            q = int(divmod(mem._ssfitness, mean)[0])
            # Account for the fact that a molecule may have been
            # yielded due to elitism.
            q -= 1 if mem in yielded else 0
            q = 0 if mem in yielded and not duplicates else q
            for x in range(q):
                yielded.add(mem)
                yield mem
                if not duplicates:
                    break

    def _stochastic_sampling_roulette(self, pop,
                                      yielded, duplicates):
        """
        Does the roulette component of stochastic sampling.

        Parameters
        ----------
        pop : :class:`list` of :class:`.MacroMolecule`
            The molecules which are being selected.

        yielded : :class:`set` of :class:`.MacroMolecule`
            Holds all previously yielded molecules.

        duplicates : :class:`bool`
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        :class:`.MacroMolecule`
            The seleceted population member.

        """

        yielded = set() if duplicates else yielded

        while True:
            valid_pop = [ind for ind in pop if ind not in yielded]

            if not valid_pop:
                break

            mean = np.mean([mem._ssfitness for mem in pop])
            rtotal = sum(divmod(ind._ssfitness, mean)[1] for
                         ind in valid_pop)
            weights = [divmod(ind._ssfitness, mean)[1] / rtotal for
                       ind in valid_pop]

            selected = np.random.choice(valid_pop, p=weights)
            yielded.add(selected) if not duplicates else None
            yield selected
