"""
Alcohol Key
===========

"""


class AlcoholKey:
    """
    Generates keys for :class:`.Alcohol` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize an :class:`.AlcoholKey` instance.

        """

        return

    def get_key(self, alcohol):
        """
        Get the key of `alcohol`.

        Parameters
        ----------
        alcohol : :class:`.Alcohol`
            The alcohol for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            alcohol.get_oxygen().get_id(),
            alcohol.get_hydrogen().get_id(),
            alcohol.get_atom().get_id(),
            tuple(sorted(alcohol.get_bonder_ids())),
            tuple(sorted(alcohol.get_deleter_ids())),
        )
