"""
Molecule Record
===============

"""


class MoleculeRecord:
    """
    An abstract base class for molecular records used by the EA.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is a default implementation provided
    purely for convenience. Subclasses can freely ignore or
    override this implementation.

    """

    def __init__(self, molecule):
        """
        Initialize a :class:`.MoleculeRecord` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule the record holds.

        """

        self._molecule = molecule
        self._fitness_value = None
        self._normalized_fitness_value = None

    def get_molecule(self):
        """
        Get the molecule held by the record.

        Returns
        -------
        :class:`.Molecule`
            The molecule held by the record.

        """

        return self._molecule

    def get_fitness_value(self, normalized=True):
        """
        Get the fitness value of the molecule in the record.

        Parameters
        ----------
        normalized : :class:`bool`, optional
            Toggles the return of the normalized vs unnormalized
            fitness value. The unnormalized fitness value is
            guaranteed to be constant for the same molecule
            across generations, while the normalized one is allowed
            to change.

        Returns
        -------
        :class:`float`
            If `normalized` is ``True`` and a fitness value has been
            assigned, it is guaranteed to be a float.

        :class:`object`
            If `normalized` is ``False`` and a fitness value has been
            assigned,any object may be returned.

        None : :class:`NoneType`
            If a fitness value has not been assigned to the record.

        """

        return (
            self._normalized_fitness_value
            if normalized
            else self._fitness_value
        )

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.MoleculeRecord`
            The clone. Has the same type as the original record.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._molecule = self._molecule
        clone._fitness_value = self._fitness_value
        self._normalized_fitness_value = self._normalized_fitness_value

    def with_fitness_value(self, fitness_value, normalized=True):
        """
        Return a clone holding a different fitness value.

        Parameters
        ----------
        fitness_value : :class:`object` or :class:`float`
            The fitness value of the clone. If `normalized` is
            ``True``, this value must be a :class:`float`.

        normalized : :class:`bool`, optional
            Toggles if the normalized or unnormalized fitness value is
            being set.

        Returns
        -------
        :class:`.MoleculeRecord`
            The clone. Has the same type as the original record.

        """

        return self.clone()._with_fitness_value(
            fitness_value=fitness_value,
            normalized=normalized,
        )

    def _with_fitness_value(self, fitness_value, normalized):
        if normalized:
            self._normalized_fitness_value = fitness_value
        else:
            self._fitness_value = fitness_value
        return self
