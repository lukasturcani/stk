import stk
import rdkit.Chem.AllChem as rdkit


def is_equivalent_atom(atom1, atom2):
    assert atom1 is not atom2
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.get_mass() == atom2.get_mass()


_periodic_table = rdkit.GetPeriodicTable()

atomic_numbers = {
    atomic_number: stk.__dict__[_periodic_table.GetElementSymbol]
    for atomic_number in range(1, 119)
}

atomic_masses = {
    atomic_number: _periodic_table.GetAtomicWeight(atomic_number)
    for atomic_number in range(1, 119)
}
