import itertools
import logging
import typing
from collections.abc import Callable, Iterator, Sequence

from stk._internal.key_makers.molecule import MoleculeKeyMaker

from ..batch import Batch
from .utilities.yielded_batches import YieldedBatches

T = typing.TypeVar("T")
logger = logging.getLogger(__name__)


class Selector(typing.Generic[T]):
    """
    An abstract base class for selectors.

    Selectors select batches of molecules from a population.
    Each batch is selected based on its fitness. The fitness of a
    batch is the sum of all fitness values of the molecules in the
    batch. Batches may be of size 1.

    Notes:

        You might notice that some of the public methods of this abstract
        base class are implemented. This is purely for convenience when
        implementing subclasses. The implemented public methods are simply
        default implementations, which can be safely ignored or overridden,
        when implementing subclasses. Any private methods are
        implementation details of these default implementations.

        *The Default Implementation*

        This section is only of use to people who want to add a new
        :class:`.Selector` subclass, and want to make use of the default
        implementation to make this job easier.

        When using the default implementation you do not need to
        implement :meth:`.select`, which is already provided, but instead
        :meth:`._select_from_batches` needs to be implemented. What the
        default implementation provides, is code, which does the batching
        of a `population` for you, which means you only have to worry
        about implementing the selection algorithm, which works on batches
        directly.

        The default implementation also automatically updates a
        :class:`.YieldedBatches` object for you, so that you can keep track
        of which batches have already been yielded, in case you want to
        prevent duplicate selection of batches or molecule records. Though
        whether you want to make use of this will depend on the nature of
        your selection algorithm.

    See Also:

        :class:`.Batch`

    Examples:

        *Subclass Implementation*

        The source code of the classes listed in :mod:`.selector` can serve
        as good examples.

    """

    def __init__(
        self,
        key_maker: MoleculeKeyMaker,
        fitness_modifier: Callable[[Sequence[T]], dict[T, float]],
        batch_size: int,
    ) -> None:
        """
        Parameters:

            key_maker:
                Used to get the keys of molecules, which are used to
                determine if two molecule records are duplicates of each
                other.

            fitness_modifier:
                Takes the `population` on which :meth:`.select` is called
                and returns a :class:`dict`, which maps records in the
                `population` to the fitness values the :class:`.Selector`
                should use.

            batch_size:
                The number of molecules yielded at once.
        """
        self._key_maker = key_maker
        self._fitness_modifier = fitness_modifier
        self._batch_size = batch_size

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ) -> Iterator[Batch]:
        """
        Yield batches of molecule records from `population`.

        Parameters:

            population : :class:`tuple` of :class:`.MoleculeRecord`
                A collection of molecules from which batches are selected.

            included_batches : :class:`set`, optional
                The identity keys of batches which are allowed to be
                yielded, if ``None`` all batches can be yielded. If not
                ``None`` only batches `included_batches` will be yielded.

            excluded_batches : class:`set`, optional
                The identity keys of batches which are not allowed to be
                yielded. If ``None``, no batch is forbidden from being
                yielded.

        Yields:

            A batch of selected molecule records.

        """

        batches = tuple(
            self._get_batches(
                population=population,
                fitness_values=self._fitness_modifier(population),
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        )

        yielded_batches = YieldedBatches(self._key_maker)
        for batch in self._select_from_batches(
            batches=batches,
            yielded_batches=yielded_batches,
        ):
            yielded_batches.update(batch)
            yield batch

        cls_name = self.__class__.__name__
        logger.debug(
            f"{cls_name} yielded {yielded_batches.get_num()} batches."
        )

    def _get_batches(
        self,
        population,
        fitness_values,
        included_batches,
        excluded_batches,
    ) -> Iterator[Batch]:
        """
        Get batches molecules from `population`.

        Parameters:

            population : :class:`tuple` of :class:`.MoleculeRecord`
                The molecule records which are to be batched.

            fitness_values : :class:`dict`
                Maps each :class:`.MoleculeRecord` in `population` to the
                fitness value which should be used for it.

            included_batches : :class:`set`, optional
                The identity keys of batches which are allowed to be
                yielded, if ``None`` all batches can be yielded. If not
                ``None`` only batches `included_batches` will be yielded.

            excluded_batches : class:`set`, optional
                The identity keys of batches which are not allowed to be
                yielded. If ``None``, no batch is forbidden from being
                yielded.

            Yields:
                A batch of molecules from `population`.

        """

        def is_included(batch):
            if included_batches is None:
                return True
            return batch.get_identity_key() in included_batches

        def is_excluded(batch):
            if excluded_batches is None:
                return False
            return batch.get_identity_key() in excluded_batches

        for records in itertools.combinations(population, self._batch_size):
            batch = Batch(
                records=records,
                fitness_values=fitness_values,
                key_maker=self._key_maker,
            )
            if is_included(batch) and not is_excluded(batch):
                yield batch

    def _select_from_batches(
        self,
        batches,
        yielded_batches,
    ) -> Iterator[Batch]:
        """
        Yield batches from `batches`.

        Parameters:
            batches : :class:`tuple` of :class:`.Batches`
                Batches from which some are selected.

            yielded_batches : :class:`.YieldedBatches`
                Keeps track of which batches have been yielded. This
                object automatically updates each time ``yield`` is called.

        Yields:
            A selected batch.

        """

        raise NotImplementedError()

    def _get_fitness_values(self, population):
        return {record: record.get_fitness_value() for record in population}
