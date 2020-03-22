"""
Boronic Acid Key Maker
======================

"""


class BoronicAcidKeyMaker:
    """
    Generates keys for :class:`.BoronicAcid` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.BoronicAcidKeyMaker` instance.

        """

        return

    def get_key(self, boronic_acid):
        """
        Get the key of `boronic_acid`.

        Parameters
        ----------
        boronic_acid : :class:`.BoronicAcid`
            The boronic acid for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'boronic_acid',
            boronic_acid.get_boron().get_id(),
            boronic_acid.get_oxygen1().get_id(),
            boronic_acid.get_hydrogen1().get_id(),
            boronic_acid.get_oxygen2().get_id(),
            boronic_acid.get_hydrogen2().get_id(),
            boronic_acid.get_atom().get_id(),
            tuple(sorted(boronic_acid.get_bonder_ids())),
            tuple(sorted(boronic_acid.get_deleter_ids())),
        )
