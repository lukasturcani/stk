"""
Alkene Key
==========

"""


class AlkeneKey:
    """
    Generates keys for :class:`.Alkene` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize an :class:`.AlkeneKey` instance.

        """

        return

    def get_key(self, alkene):
        """
        Get the key of `alkene`.

        Parameters
        ----------
        alkene : :class:`.Alkene`
            The alkene for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'alkene',
            alkene.get_carbon1().get_id(),
            alkene.get_atom1().get_id(),
            alkene.get_atom2().get_id(),
            alkene.get_carbon2().get_id(),
            alkene.get_atom3().get_id(),
            alkene.get_atom4().get_id(),
            tuple(sorted(alkene.get_bonder_ids())),
            tuple(sorted(alkene.get_deleter_ids())),
        )
