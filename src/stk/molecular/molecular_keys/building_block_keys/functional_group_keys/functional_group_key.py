"""
Functional Group Key
====================

"""


class FunctionalGroupKey:
    """
    Generates keys for :class:`.FunctionalGroup` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.FunctionalGroupKey` instance.

        """

        return

    def get_key(self, functional_group):
        """
        Get the key of `functional_group`.

        Parameters
        ----------
        functional_group : :class:`.FunctionalGroup`
            The functionalGroup for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'functional_group',
            tuple(sorted(functional_group.get_atom_ids())),
            tuple(sorted(functional_group.get_placer_ids())),
            tuple(sorted(functional_group.get_core_atom_ids())),
        )
