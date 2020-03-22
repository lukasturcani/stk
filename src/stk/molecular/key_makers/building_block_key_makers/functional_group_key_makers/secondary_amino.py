"""
Secondary Amino Key Maker
=========================

"""


class SecondaryAminoKeyMaker:
    """
    Generates keys for :class:`.SecondaryAmino` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.SecondaryAminoKey` instance.

        """

        return

    def get_key(self, secondary_amino):
        """
        Get the key of `secondary_amino`.

        Parameters
        ----------
        secondary_amino : :class:`.SecondaryAmino`
            The secondary amino functional group for which a key is
            needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'secondary_amino',
            secondary_amino.get_nitrogen().get_id(),
            secondary_amino.get_hydrogen().get_id(),
            secondary_amino.get_atom1().get_id(),
            secondary_amino.get_atom2().get_id(),
            tuple(sorted(secondary_amino.get_bonder_ids())),
            tuple(sorted(secondary_amino.get_deleter_ids())),
        )
