"""
Thioacid Key Maker
==================

"""


class ThioacidKeyMaker:
    """
    Generates keys for :class:`.Thioacid` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.ThioacidKeyMaker` instance.

        """

        return

    def get_key(self, thioacid):
        """
        Get the key of `thioacid`.

        Parameters
        ----------
        thioacid : :class:`.Thioacid`
            The thioacid for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'thioacid',
            thioacid.get_carbon().get_id(),
            thioacid.get_oxygen().get_id(),
            thioacid.get_sulfur().get_id(),
            thioacid.get_hydrogen().get_id(),
            thioacid.get_atom().get_id(),
            tuple(sorted(thioacid.get_bonder_ids())),
            tuple(sorted(thioacid.get_deleter_ids())),
        )
