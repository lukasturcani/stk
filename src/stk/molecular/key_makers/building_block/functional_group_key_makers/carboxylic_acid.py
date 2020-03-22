"""
Carboxylic Acid Key Maker
=========================

"""


class CarboxylicAcidKeyMaker:
    """
    Generates keys for :class:`.CarboxylicAcid` instances.

    """

    # Keep the empty __init__() method for the docstring.
    def __init__(self):
        """
        Initialize a :class:`.CarboxylicAcidKeyMaker` instance.

        """

        return

    def get_key(self, carboxylic_acid):
        """
        Get the key of `carboxylic_acid`.

        Parameters
        ----------
        carboxylic_acid : :class:`.CarboxylicAcid`
            The carboxylic acid for which a key is needed.

        Returns
        -------
        :class:`object`
            The key.

        """

        return (
            'carboxylic_acid',
            carboxylic_acid.get_carbon().get_id(),
            carboxylic_acid.get_oxygen1().get_id(),
            carboxylic_acid.get_oxygen2().get_id(),
            carboxylic_acid.get_hydrogen().get_id(),
            carboxylic_acid.get_atom().get_id(),
            tuple(sorted(carboxylic_acid.get_bonder_ids())),
            tuple(sorted(carboxylic_acid.get_deleter_ids())),
        )
