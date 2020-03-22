"""
Fluoro Key
==========

"""


class FluoroKey:
    """
    Generates keys for :class:`.FluoroKey` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.FluoroKey` instance.

        """

        return

    def get_key(self, fluoro):
        """
        Get the key of `fluoro`.

        Parameters
        ----------
        fluoro : :class:`.Fluoro`
            The fluoro for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'fluoro',
            fluoro.get_fluorine().get_id(),
            fluoro.get_atom().get_id(),
            tuple(sorted(fluoro.get_bonder_ids())),
            tuple(sorted(fluoro.get_deleter_ids())),
        )
