"""
Generic Functional Group Key Maker
==================================

"""


class GenericFunctionalGroupKeyMaker:
    """
    Generates keys for :class:`.GenericFunctionalGroup` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.GenericFunctionalGroupKeyMaker` instance.

        """

        return

    def get_key(self, generic_functional_group):
        """
        Get the key of `generic_functional_group`.

        Parameters
        ----------
        generic_functional_group : :class:`.GenericFunctionalGroup`
            The generic_functional_group for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'generic_functional_group',
            tuple(sorted(generic_functional_group.get_atom_ids())),
            tuple(sorted(generic_functional_group.get_bonder_ids())),
            tuple(sorted(generic_functional_group.get_deleter_ids())),
        )
