class FunctionalGroup:
    """
    An abstract base class for functional groups.

    """

    def get_atoms(self):
        raise NotImplementedError()

    def get_atom_ids(self):
        raise NotImplementedError()

    def get_bonders(self):
        raise NotImplementedError()

    def get_bonder_ids(self):
        raise NotImplementedError()

    def get_deleters(self):
        raise NotImplementedError()

    def get_deleter_ids(self):
        raise NotImplementedError()
