"""
Iodo Key Maker
==============

"""


class IodoKeyMaker:
    """
    Generates keys for :class:`.Iodo` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.IodoKeyMaker` instance.

        """

        return

    def get_key(self, iodo):
        """
        Get the key of `iodo`.

        Parameters
        ----------
        iodo : :class:`.Iodo`
            The iodo for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'iodo',
            iodo.get_iodine().get_id(),
            iodo.get_atom().get_id(),
            tuple(sorted(iodo.get_bonder_ids())),
            tuple(sorted(iodo.get_deleter_ids())),
        )
