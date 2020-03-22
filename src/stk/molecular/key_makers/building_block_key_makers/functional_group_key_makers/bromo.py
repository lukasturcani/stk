"""
Bromo Key Maker
===============

"""


class BromoKeyMaker:
    """
    Generates keys for :class:`.Bromo` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.BromoKeyMaker` instance.

        """

        return

    def get_key(self, bromo):
        """
        Get the key of `bromo`.

        Parameters
        ----------
        bromo : :class:`.Bromo`
            The bromo group for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'bromo',
            bromo.get_bromine().get_id(),
            bromo.get_atom().get_id(),
            tuple(sorted(bromo.get_bonder_ids())),
            tuple(sorted(bromo.get_deleter_ids())),
        )
