"""
Selector
========

#. :class:`.AboveAverage`
#. :class:`.Best`
#. :class:`.FilterBatches`
#. :class:`.FilterMoleculeRecords`
#. :class:`.RemoveBatches`
#. :class:`.RemoveMoleculeRecords`
#. :class:`.Roulette`
#. :class:`.StochasticUniversalSampling`
#. :class:`.Tournament`
#. :class:`.Worst`

"""

import logging
import itertools as it

from ..batch import Batch

logger = logging.getLogger(__name__)


class Selector:
    """
    An abstract base class for selectors.

    Selectors select batches of molecules from a population.
    Each batch is selected based on its fitness. The fitness of a
    batch is the sum of all fitness values of the molecules in the
    batch. Batches may be of size 1.

    See Also
    --------
    :class:`.Batch`

    Examples
    --------
    *Subclass Implementation*

    The source code of the classes listed in :mod:`.selector` can serve
    as good examples.

    """

    def __init__(self, key_maker, fitness_modifier):
        self._key_maker = key_maker
        self._fitness_modifier = fitness_modifier

    def select(
        self,
        population,
        included_batches=None,
        excluded_batches=None,
    ):
        """
        Select batches of molecules from `population`.

        Parameters
        ----------
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

        Yields
        ------
        :class:`Batch` of :class:`.MoleculeRecord`
            A batch of selected molecule records.

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

        Yields
        ------
        :class:`.Batch`
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

        for records in it.combinations(population, self._batch_size):
            batch = Batch(
                records=records,
                fitness_values=fitness_values,
                key_maker=self._key_maker,
            )
            if is_included(batch) and not is_excluded(batch):
                yield batch

    def _select_from_batches(self, batches, yielded):
        """

        """

        raise NotImplementedError()
