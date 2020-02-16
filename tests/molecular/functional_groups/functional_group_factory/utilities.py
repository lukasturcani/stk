import itertools as it


class _TestCase:
    """
    A test case.

    Attributes
    ----------
    factory : :class:`.FunctionalGroupFactory`
        The factory being tested.

    molecule : :class:`.Molecule`
        The molecule passed to
        :meth:`.FunctionalGroupFactory.get_functional_groups`.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The expected functional groups.

    """

    def __init__(self, factory, molecule, functional_groups):
        self.factory = factory
        self.molecule = molecule
        self.functional_groups = functional_groups


def atom_id(atom):
    return atom.get_id()


def are_same_id_sequences(ids1, ids2):
    for id1, id2 in it.zip_longest(ids1, ids2):
        assert id1 == id2


def are_clone_sequences(atoms1, atoms2):
    """
    Test if `atoms1` and `atoms2` are clones of each other.

    """

    for a1, a2 in it.zip_longest(atoms1, atoms2):
        assert a1 is not a2
        assert a1.get_id() == a2.get_id()
        assert a1.get_charge() == a2.get_charge()
        assert a1.__class__ is a2.__class__
