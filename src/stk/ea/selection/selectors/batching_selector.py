import itertools as it
import logging

from ..batch import Batch

logger = logging.getLogger(__name__)


class Batcher:
    """
    Creates batches.

    """

    @staticmethod
    def _return_fitness_values(population):
        return population.get_fitness_values()



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
