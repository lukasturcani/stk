"""
Thiol Key
=========

"""


class ThiolKey:
    """
    Generates keys for :class:`.ThiolKey` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.ThiolKey` instance.

        """

        return

    def get_key(self, thiol):
        """
        Get the key of `thiol`.

        Parameters
        ----------
        thiol : :class:`.Thiol`
            The thiol for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'thiol',
            thiol.get_sulfur().get_id(),
            thiol.get_hydrogen().get_id(),
            thiol.get_atom().get_id(),
            tuple(sorted(thiol.get_bonder_ids())),
            tuple(sorted(thiol.get_deleter_ids())),
        )
