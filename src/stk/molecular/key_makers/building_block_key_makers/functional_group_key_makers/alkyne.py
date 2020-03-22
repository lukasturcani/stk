"""
Alkyne Key
==========

"""


class AlkyneKey:
    """
    Generates keys for :class:`.Alkyne` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize an :class:`.AlkyneKey` instance.

        """

    def get_key(self, alkyne):
        """
        Get the key of `alkyne`.

        Parmaeters
        ----------
        alkyne : :class:`.Alkyne`
            The alkyne for which is a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'alkyne',
            alkyne.get_carbon1().get_id(),
            alkyne.get_atom1().get_id(),
            alkyne.get_carbon2().get_id(),
            alkyne.get_atom2().get_id(),
            tuple(sorted(alkyne.get_bonder_ids())),
            tuple(sorted(alkyne.get_deleter_ids())),
        )
