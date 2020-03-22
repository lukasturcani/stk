"""
Difluoro Key Maker
==================

"""


class DifluoroKeyMaker:
    """
    Generates keys for :class:`.Difluoro` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.DifluoroKeyMaker` instance.

        """

        return

    def get_key(self, difluoro):
        """
        Get the key of `difluoro`.

        Parameters
        ----------
        difluoro : :class:`.Difluoro`
            The difluoro for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'difluoro',
            difluoro.get_fluorine1().get_id(),
            difluoro.get_atom1().get_id(),
            difluoro.get_fluorine2().get_id(),
            difluoro.get_atom2().get_id(),
            tuple(sorted(difluoro.get_bonder_ids())),
            tuple(sorted(difluoro.get_deleter_ids())),
        )
