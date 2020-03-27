
class ShiftUp(_FilteringNormalizer, FitnessNormalizer):
    """
    Shifts negative fitness values to be positive.

    Assume you have a vector-valued fitness value, where each number
    represents a different property of the molecule

    .. code-block:: python

        {mol1: [1, -10, 1]}

    One way to convert the vector-valued fitness value into a
    scalar fitness value is by summing the elements, and the result in
    this case would be ``-8``. Clearly this doesn't work, because the
    resulting fitness value is not a positive number. To fix this,
    the ``-10`` should be shifted to a positive value.

    :class:`.ShiftUp` finds the minimum value of each element in the
    vector-valued fitness value across the entire population, and for
    element where this minimum value is less than ``0``, shifts up
    the element value for every molecule in the population, so that the
    minimum value in the entire population is ``1``.

    For example, take a population with the vector-valued fitness
    values

    .. code-block:: python

        fitness_values = {
            mol1: [1, -5, 5],
            mol2: [3, -10, 2],
            mol3: [2, 20, 1],
        }

    After normalization the fitness values will be.

    .. code-block:: python

        normalized  = {
            mol1: [1, 6, 5],
            mol2: [3, 1, 2],
            mol3: [2, 31, 1],
        }

    This :class:`.ShiftUp` also works when the fitness value is a
    single value.

    Examples
    --------
    .. code-block:: python

        # Create the normalizer.
        shifter = ShiftUp()

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, -2, 3], mol2: [4, 5, -6], mol3: [7, 8, 9]}.
        normalized = shifter.normalize(pop)

        # normalized is
        # {mol1: [1, 1, 10], mol2: [4, 8, 1], mol3: [7, 11, 16]}.

    """

    def __init__(self, filter=lambda population, mol: True):
        """
        Initialize a :class:`.ShiftUp` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a two parameters, first is a :class:`.EAPopulation`
            and the second is a :class:`.Molecule`, and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.
            The :class:`.EAPopulation` on which :meth:`normalize` is
            called is passed as the first argument while the second
            argument will be passed every :class:`.Molecule` in it.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        # filtered is iterated through multiple times.
        filtered = list(filtered)

        # Get all the fitness arrays in a matrix.
        fmat = np.array([fitness_values[mol] for mol in filtered])

        # Get the minimum value of each element across the population.
        # keepdims ensures that np.min returns a 1-D array, because
        # it will be True if fitness values are scalar and False if
        # they are array-valued.
        mins = np.min(fmat, axis=0, keepdims=len(fmat.shape) == 1)

        # Convert all elements in mins which are not to be shifted to 0
        # and make the shift equal to the minimum value + 1.
        shift = np.zeros(len(mins))
        for i, min_ in enumerate(mins):
            if min_ <= 0:
                shift[i] = 1 - min_

        for mol in filtered:
            yield mol, fitness_values[mol] + shift


