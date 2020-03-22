"""
Diol Key
========

"""


class DiolKey:
    """
    Generates keys for :class:`.Diol` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.DiolKey` instance.

        """

        return

    def get_key(self, diol):
        """
        Get the key of `diol`.

        Parameters
        ----------
        diol : :class:`.Diol`
            The diol for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'diol',
            diol.get_atom1().get_id(),
            diol.get_oxygen1().get_id(),
            diol.get_hydrogen1().get_id(),
            diol.get_atom2().get_id(),
            diol.get_oxygen2().get_id(),
            diol.get_hydrogen2().get_id(),
            tuple(sorted(diol.get_bonder_ids())),
            tuple(sorted(diol.get_deleter_ids())),
        )
