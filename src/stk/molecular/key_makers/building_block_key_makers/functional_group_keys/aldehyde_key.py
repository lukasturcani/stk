"""
Aldehyde Key
============

"""


class AldehydeKey:
    """
    Generates keys for :class:`.Aldehyde` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.AldehydeKey` instance.

        """

        return

    def get_key(self, aldehyde):
        """
        Get the key of `aldehyde`.

        Parameters
        ----------
        aldehyde : :class:`.Aldehyde`
            The aldehyde for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'aldehyde',
            aldehyde.get_carbon().get_id(),
            aldehyde.get_oxygen().get_id(),
            aldehyde.get_hydrogen().get_id(),
            aldehyde.get_atom().get_id(),
            tuple(sorted(aldehyde.get_bonder_ids())),
            tuple(sorted(aldehyde.get_deleter_ids())),
        )
