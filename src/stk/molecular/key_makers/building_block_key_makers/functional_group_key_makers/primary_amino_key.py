"""
Primary Amino Key
=================

"""


class PrimaryAminoKey:
    """
    Generates keys for :class:`.PrimaryAmino` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.PrimaryAminoKey` instance.

        """

        return

    def get_key(self, primary_amino):
        """
        Get the key of `primary_amino`.

        Parameters
        ----------
        primary_amino : :class:`.PrimaryAmino`
            The primary amino functional group for which a key is
            needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'primary_amino',
            primary_amino.get_nitrogen().get_id(),
            primary_amino.get_hydrogen1().get_id(),
            primary_amino.get_hydrogen2().get_id(),
            primary_amino.get_atom().get_id(),
            tuple(sorted(primary_amino.get_bonder_ids())),
            tuple(sorted(primary_amino.get_deleter_ids())),
        )
