#!/usr/bin/env python3
import os
import numpy as np
from copy import deepcopy

from .utilities import (
    discrete_molecules, decipher_atom_key, molecular_weight, center_of_mass,
    max_dim, void_diameter, opt_void_diameter, void_volume, shift_com,
    find_windows, create_supercell, make_JSON_serializable
)
from .io_tools import Input, Output


class _NotAModularSystem(Exception):
    def __init__(self, message):
        self.message = message


class Molecule(object):
    def __init__(self, mol, system_name, mol_id):
        self._Output = Output()
        self.mol = mol
        self.no_of_atoms = len(mol['elements'])
        self.elements = mol['elements']
        if 'atom_ids' in mol.keys():
            self.atom_ids = mol['atom_ids']
        self.coordinates = mol['coordinates']
        self.parent_system = system_name
        self.molecule_id = mol_id
        self.properties = {}

    @classmethod
    def load_rdkit_mol(cls, mol, system_name='rdkit', mol_id=0):
        return cls(Input().load_rdkit_mol(mol), system_name, mol_id)

    def full_analysis(self, ncpus=1, **kwargs):
        results = {
            'no_of_atoms': self.no_of_atoms,
            'mol_weight': self.molecular_weight,
            'COM': self.calculate_COM(),
            'maximum_diameter': self.calculate_maximum_diameter(),
            'void_diameter': self.calculate_void_diameter(),
            'void_volume': self.calculate_void_volume(),
            'void_diameter_opt': self.calculate_void_diameter_opt(**kwargs),
            'void_volume_opt': self.calculate_void_volume_opt(**kwargs),
            'windows': self.calculate_windows(
                ncpus=ncpus, **kwargs),
        }
        return (results)

    def calculate_COM(self):
        self.centre_of_mass = center_of_mass(self.elements, self.coordinates)
        self.properties['centre_of_mass'] = self.centre_of_mass
        return self.centre_of_mass

    def calculate_maximum_diameter(self):
        self.maxd_atom_1, self.maxd_atom_2, self.maximum_diameter = max_dim(
            self.elements, self.coordinates)
        self.properties['maximum_diameter'] = {
            'diameter': self.maximum_diameter,
            'atom_id_1': int(self.maxd_atom_1),
            'atom_id_2': int(self.maxd_atom_2),
        }
        return self.maximum_diameter

    def calculate_void_diameter(self):
        self.void_diameter, self.void_closest_atom = void_diameter(
            self.elements, self.coordinates)
        self.properties['void_diameter'] = {
            'void_diameter': self.void_diameter,
            'atom_id_1': int(self.void_closest_atom),
        }
        return self.void_diameter

    def calculate_void_volume(self):
        self.void_volume = void_volume(self.calculate_void_diameter() / 2)
        self.properties['void_volume'] = self.void_volume
        return self.void_volume

    def calculate_void_diameter_opt(self, **kwargs):
        (self.void_diameter_opt, self.void_opt_closest_atom,
         self.void_opt_COM) = opt_void_diameter(self.elements,
                                                self.coordinates, **kwargs)
        self.properties['void_diameter_opt'] = {
            'void_diameter': self.void_diameter_opt,
            'atom_id_1': int(self.void_opt_closest_atom),
            'void_COM': self.void_opt_COM,
        }
        return self.void_diameter_opt

    def calculate_void_volume_opt(self, **kwargs):
        self.void_volume_opt = void_volume(
            self.calculate_void_diameter_opt(**kwargs) / 2)
        self.properties['void_void_volume_opt'] = self.void_volume_opt
        return self.void_volume_opt

    def calculate_windows(self, **kwargs):
        windows = find_windows(self.elements, self.coordinates, **kwargs)
        if 'output' in kwargs:
            if kwargs['output'] == 'windows':
                self.properties['windows'] = {'windows_diameters': windows, }
        else:
            self.properties['windows'] = {
                'windows_diameters': windows[0],
                'windows_coms': windows[1],
            }
        return windows

    def shift_to_origin(self, **kwargs):
        self.coordinates = shift_com(self.elements, self.coordinates, **kwargs)
        self._update()

    def molecular_weight(self, mol):
        self.molecular_weight = molecular_weight(mol['elements'])
        self.properties['windows'] = self.molecular_weight
        return self.molecular_weight

    def save_molecule_json(self, filepath=None, molecular=False, **kwargs):
        # We pass a copy of the properties dictionary.
        dict_obj = deepcopy(self.properties)
        # If molecular data is also required we update the dictionary.
        if molecular is True:
            dict_obj.update(self.mol)
        # We make sure it is JSON serializable.
        dict_obj = make_JSON_serializable(dict_obj)
        # If no filepath is provided we create one.
        if filepath is None:
            filepath = "_".join(
                (str(self.parent_system), str(self.molecule_id))
            )
            filepath = '/'.join((os.getcwd(), filepath))
        # Dump the dictionary to json file.
        self._Output.dump2json(dict_obj, filepath, **kwargs)

    def save_molecule(self, filepath=None, **kwargs):
        # If no filepath is provided we create one.
        if filepath is None:
            filepath = "_".join(
                (str(self.parent_system), str(self.molecule_id))
            )
            filepath = '/'.join((os.getcwd(), filepath))
            filepath = '.'.join((filepath, 'pdb'))
        # Check if there is an 'atom_ids' keyword in the self.mol dict.
        # Otherwise pass to the dump2file atom_ids='elements'.
        if 'atom_ids' not in self.mol.keys():
            atom_ids = 'elements'
        else:
            atom_ids = 'atom_ids'
        # Dump molecule into a file.
        self._Output.dump2file(self.mol, filepath, atom_ids=atom_ids, **kwargs)

    def _update(self):
        self.mol['coordinates'] = self.coordinates
        self.calculate_COM()
        self.calculate_void_diameter_opt()


class MolecularSystem(object):
    def __init__(self):
        self._Input = Input()
        self._Output = Output()
        self.system_id = 0

    @classmethod
    def load_file(cls, file_path):
        obj = cls()
        obj.system = obj._Input.load_file(file_path)
        obj.filename = os.path.basename(file_path)
        obj.system_id = obj.filename.split(".")[0]
        obj.name, ext = os.path.splitext(obj.filename)
        return obj

    @classmethod
    def load_rdkit_mol(cls, mol):
        obj = cls()
        obj.system = obj._Input.load_rdkit_mol(mol)
        return obj

    @classmethod
    def load_system(cls, dict_):
        obj = cls()
        obj.system = dict_
        return obj

    def reconstruct_system(self, **kwargs):
        # First we create a 3x3x3 supercell with the initial unit cell in the
        # centre and the 26 unit cell translations around to provide all the
        # atom positions necessary for the molecules passing through periodic
        # boundary reconstruction step.
        supercell = create_supercell(self.system, **kwargs)
        discrete = discrete_molecules(supercell)
        return discrete

    def swap_atom_keys(self, swap_dict, dict_key='atom_ids'):
        """
        Swap atom_key for atom_key in system's 'elements' array.

        Parameters
        ----------
        swap_dict: dict
            A dictionary containg atom keys (dictionary's keys) to be swapped
            with corresponding atom keys (dictionary's keys arguments).

        dict_key: str (default='elements')
            A key in MolecularSystem().system dictionary to perform the
            atom keys swapping operation on.

        Modifies
        --------
        system['elements']: array
            Replaces every occurance of dictionary key with dictionary key's
            argument in the 'elements' array of the MolecularSystem().system's
            dictionary.

        Returns
        -------
        None: NoneType

        """
        # Similar situation to the one from decipher_atom_keys function.
        if 'atom_ids' not in self.system.keys():
            dict_key = 'elements'
        for atom_key in range(len(self.system[dict_key])):
            for key in swap_dict.keys():
                if self.system[dict_key][atom_key] == key:
                    self.system[dict_key][atom_key] = swap_dict[key]

    def decipher_atom_keys(self, forcefield='DLF', dict_key='atom_ids'):
        """
        Decipheres forcefield's keys for their periodic elements equivalents.

        This function runs decipher_atom_key() function on every atom key in
        the system['atom_keys'] array and substitutes the return of this
        function at the corresponding position of system['elements'] array.

        The supported forcfields is OPLS and also the DL_F notation
        (see User's Guide) with keywords allowed:
        'OPLS', 'OPLS2005', 'OPLSAA', 'OPLS3' and 'DLF', 'DL_F'.

        Parameters
        ----------
        forcefield: str
            The forcefield used to decipher the atom keys. This parameter is
            not case sensitive.

        Modifies
        --------
        system['elements']
            It substitutes the string objects in this array for the return
            string of the decipher_atom_key() for each atom key in
            system['atom_keys'] array equvalent.

        Returns
        -------
        None: NoneType

        """
        # In case there is no 'atom_ids' key we try 'elements'. This is for
        # XYZ and MOL files mostly. But, we keep the dict_key keyword for
        # someone who would want to decipher 'elements' even if 'atom_ids' key
        # is present in the system's dictionary.
        if 'atom_ids' not in self.system.keys():
            dict_key = 'elements'
        # I do it on temporary object so that it only finishes when successful
        temp = deepcopy(self.system[dict_key])
        for element in range(len(temp)):
            temp[element] = "{0}".format(
                decipher_atom_key(
                    temp[element], forcefield=forcefield))
        self.system['elements'] = temp

    def make_modular(self, supercell=False):
        if supercell is True:
            supercell_333 = create_supercell(self.system)
        else:
            supercell_333 = None
        dis = discrete_molecules(self.system, supercell=supercell_333)
        self.no_of_discrete_molecules = len(dis)
        # elements_ = []
        # coordinates_ = []
        self.molecules = {}
        for i in range(len(dis)):
            self.molecules[i] = Molecule(dis[i], self.name, i)
            # [elements_.append(i) for i in dis[i]['elements']]
            # [coordinates_.append(list(i)) for i in dis[i]['coordinates']]
        # self.system['modular_elements'] = np.array(elements_)
        # self.system['modular_coordinates'] = np.array(coordinates_)

    def system_to_molecule(self):
        return Molecule(self.system, self.system_id, 0)

    def save_system(self, filepath=None, modular=False, **kwargs):
        # If no filepath is provided we create one.
        if filepath is None:
            filepath = '/'.join((os.getcwd(), str(self.system_id)))
            filepath = '.'.join((filepath, 'pdb'))
        # If modular is True substitute the molecular data for modular one.
        system_dict = deepcopy(self.system)
        if modular is True:
            elements = np.array([])
            atom_ids = np.array([])
            coor = np.array([]).reshape(0, 3)
            for mol_ in self.molecules:
                mol = self.molecules[mol_]
                elements = np.concatenate((elements, mol.mol['elements']))
                atom_ids = np.concatenate((atom_ids, mol.mol['atom_ids']))
                coor = np.concatenate((coor, mol.mol['coordinates']), axis=0)
            system_dict['elements'] = elements
            system_dict['atom_ids'] = atom_ids
            system_dict['coordinates'] = coor
        # Check if there is an 'atom_ids' keyword in the self.mol dict.
        # Otherwise pass to the dump2file atom_ids='elements'.
        # This is mostly for XYZ files and not deciphered trajectories.
        if 'atom_ids' not in system_dict.keys():
            atom_ids = 'elements'
        else:
            atom_ids = 'atom_ids'
        # Dump system into a file.
        self._Output.dump2file(
            system_dict, filepath, atom_ids=atom_ids, **kwargs)

    def save_system_json(self, filepath=None, modular=False, **kwargs):
        # We pass a copy of the properties dictionary.
        dict_obj = deepcopy(self.system)
        # We make sure it is JSON serializable.
        dict_obj = make_JSON_serializable(dict_obj)
        # In case we want a modular system.
        if modular is True:
            try:
                if self.molecules:
                    pass
            except AttributeError:
                raise _NotAModularSystem(
                    "This system is not modular. Please, run first the "
                    "make_modular() function of this class.")
            dict_obj = {}
            for molecule in self.molecules:
                mol_ = self.molecules[molecule]
                mol_dict = make_JSON_serializable(mol_.mol)
                dict_obj[molecule] = mol_dict
        # If no filepath is provided we create one.
        if filepath is None:
            filepath = '/'.join((os.getcwd(), str(self.system_id)))
        # Dump the dictionary to json file.
        self._Output.dump2json(dict_obj, filepath, **kwargs)
