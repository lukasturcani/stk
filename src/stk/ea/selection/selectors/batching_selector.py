from .selector import Selector


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


