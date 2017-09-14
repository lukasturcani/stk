#!/usr/bin/env python3
"""
Module defining classes and functions for input/output processing.

TO-DO-LIST:
-----------
-Add appropriate comments and descriptions to all functions
"""

import os
import json
import numpy as np

from .utilities import decipher_atom_key, unit_cell_to_lattice_array


class _CorruptedPDBFile(Exception):
    def __init__(self, message):
        self.message = message


class _CorruptedXYZFile(Exception):
    def __init__(self, message):
        self.message = message


class _FileAlreadyExists(Exception):
    def __init__(self, message):
        self.message = message


class _NotADictionary(Exception):
    def __init__(self, message):
        self.message = message


class _FileTypeError(Exception):
    def __init__(self, message):
        self.message = message


class Input(object):
    """ Class for loading and processing input files. """

    def __init__(self):
        self._load_funcs = {
            '.xyz': self._read_xyz,
            '.pdb': self._read_pdb,
            '.mol': self._read_mol,
        }

    def load_file(self, filepath):
        """
        This function opens any type of a readable file and decompose
        the file object into a list, for each line, of lists containing
        splitted line strings using space as a spacer.

        Parameters
        ----------
        filepath : str
            The full path or a relative path to any type of file.

        Returns
        -------
        dict()
            Returns a dictionary containing the molecular information
            extracted from the input files. This information will
            vary with file type and information stored in it.
            The data is sorted into lists that contain one feature
            for example key atom_id: [atom_id_1, atom_id_2]
            Over the process of analysis this dictionary will be updated
            with new data.
        """
        self.file_path = filepath
        _, self.file_type = os.path.splitext(filepath)
        _, self.file_name = os.path.split(filepath)
        with open(filepath) as ffile:
            self.file_content = ffile.readlines()

        return (self._load_funcs[self.file_type]())

    def load_rdkit_mol(self, mol):
        """
        Return molecular data from mol object.

        Parameters
        ----------
        mol : Mol
            A molecule object from RDKit.

        Returns
        -------
        dictionary
            A dictionary with ``elements`` and ``coordinates`` as keys
            containing molecular data extracted from RDKit Mol object.

        """
        self.system = {
            'elements': np.empty(
                mol.GetNumAtoms(), dtype=str),
            'coordinates': np.empty((mol.GetNumAtoms(), 3))
        }
        for atom in mol.GetAtoms():
            atom_id = atom.GetIdx()
            atom_sym = atom.GetSymbol()
            x, y, z = mol.GetConformer().GetAtomPosition(atom_id)
            self.system['elements'][atom_id] = atom_sym
            self.system['coordinates'][atom_id] = x, y, z
        return self.system

    def _read_xyz(self):
        """"""
        try:
            self.system = dict()
            self.file_remarks = self.file_content[1]
            self.system['elements'] = np.array(
                [i.split()[0] for i in self.file_content[2:]])
            self.system['coordinates'] = np.array(
                [[float(j[0]), float(j[1]), float(j[2])]
                 for j in [i.split()[1:] for i in self.file_content[2:]]])
            return self.system
        except IndexError:
            raise _CorruptedXYZFile(
                "The XYZ input file is either corrupted in some way"
                " empty line at the end etc.) or if it is an XYZ trajectory"
                " file, please use Trajectory class.")

    def _read_pdb(self):
        """"""
        if sum([i.count('END ') for i in self.file_content]) > 1:
            raise _CorruptedPDBFile(
                "Multiple 'END' keywords were found in the PDB file."
                "If this is a trajectory, use Trajectory class, or fix it.")
        self.system = dict()
        self.system['remarks'] = [
            i for i in self.file_content if i[:6] == 'REMARK'
        ]
        self.system['unit_cell'] = np.array([
            float(x)
            for i in self.file_content for x in
            [i[6:15], i[15:24], i[24:33], i[33:40], i[40:47], i[47:54]]
            if i[:6] == 'CRYST1'
        ])
        if self.system['unit_cell'].any():
            self.system['lattice'] = unit_cell_to_lattice_array(self.system[
                'unit_cell'])
        self.system['atom_ids'] = np.array(
            [
                i[12:16].strip() for i in self.file_content
                if i[:6] == 'HETATM' or i[:6] == 'ATOM  '
            ],
            dtype='<U8')
        self.system['elements'] = np.array(
            [
                i[76:78].strip() for i in self.file_content
                if i[:6] == 'HETATM' or i[:6] == 'ATOM  '
            ],
            dtype='<U8')
        self.system['coordinates'] = np.array(
            [[float(i[30:38]), float(i[38:46]), float(i[46:54])]
             for i in self.file_content
             if i[:6] == 'HETATM' or i[:6] == 'ATOM  '])
        return self.system

    def _read_mol(self):
        """"""
        self.system = dict()
        if self.file_content[2] != '\n':
            self.system['remarks'] = self.file_content[2]
        file_body = [i.split() for i in self.file_content]
        elements = []
        coordinates = []
        atom_data = False
        for line in file_body:
            if len(line) > 2:
                if line[2] == 'END' and line[3] == 'ATOM':
                    atom_data = False
                if atom_data is True:
                    elements.append(line[3])
                    coordinates.append(line[4:7])
                if line[2] == 'BEGIN' and line[3] == 'ATOM':
                    atom_data = True
        self.system['elements'] = np.array(elements)
        self.system['coordinates'] = np.array(coordinates, dtype=float)
        return self.system


class Output(object):
    """ Class for saving molecular structures and/or properties to files. """

    def __init__(self):
        self.cwd = os.getcwd()
        self._save_funcs = {
            'xyz': self._save_xyz,
            'pdb': self._save_pdb,
        }

    def dump2json(self, obj, filepath, override=False, **kwargs):
        """
        This function dumps a dictionary into a json file.

        It uses the json.dump() function. All kwargs will also be passed to
        this function.
        """
        # We make sure that the object passed by the user is a dictionary.
        if isinstance(obj, dict):
            pass
        else:
            raise _NotADictionary(
                "The function only accepts object of type: dictionary.")
        # We check if the filepath has a json extenstion, if not we add it.
        if str(filepath[-4:]) == 'json':
            pass
        else:
            filepath = ".".join((str(filepath), "json"))
        # First we check if the file already exists. If yes and the override
        # keyword is False (default), we will raise an exception. Otherwise
        # the file will be overwritten.
        if override is False:
            if os.path.isfile(filepath):
                raise _FileAlreadyExists(
                    "The file you want to dump the json into alreasy exists."
                    " Change the filepath, or set 'override' to True.")
        # We dump the object to the json file. Additional kwargs can be passed.
        with open(filepath, 'w+') as json_file:
            json.dump(obj, json_file, **kwargs)

    def dump2file(self, obj, filepath, override=False, **kwargs):
        # First we check if the file already exists. If yes and the override
        # keyword is False (default), we will raise an exception. Otherwise
        # the file will be overwritten.
        if override is False:
            if os.path.isfile(filepath):
                raise _FileAlreadyExists(
                    "The file {0} alreasy exists. Change the filepath, "
                    "or set 'override' to True.".format(filepath))
        if str(filepath[-3:]) not in self._save_funcs.keys():
            raise _FileTypeError(
                "This file type is not supported for saving molecular systems"
                " and molecules. Please use XYZ or PDB."
                " Offending file type: {0}".format(str(filepath[-3:])))
        self._save_funcs[str(filepath[-3:])](obj, filepath, **kwargs)

    def _save_xyz(self, system, filepath, **kwargs):
        """"""
        # Initial settings.
        settings = {
            'elements': 'elements',
            'coordinates': 'coordinates',
            'remark': " ",
            'decipher': False,
            'forcefield': None,
        }
        settings.update(kwargs)
        # Extract neccessary data.
        elements = system['elements']
        coordinates = system['coordinates']
        if settings['decipher'] is True:
            elements = np.array([
                decipher_atom_key(
                    key, forcefield=settings['forcefield']) for key in elements
            ])
        string = '{0:0d}\n{1}\n'.format(len(elements), str(settings['remark']))
        for i, j in zip(elements, coordinates):
            string += '{0} {1:.2f} {2:.2f} {3:.2f}\n'.format(i, *j)
        with open(filepath, 'w+') as file_:
            file_.write(string)

    def _save_pdb(self, system, filepath, **kwargs):
        """"""
        settings = {
            'atom_ids': 'atom_ids',
            'elements': 'elements',
            'coordinates': 'coordinates',
            'cryst': 'unit_cell',
            'connect': None,
            'remarks': None,
            'space_group': None,
            'resName': "MOL",
            'chainID': 'A',
            'resSeq': 1,
            'decipher': False,
            'forcefield': None,
        }
        settings.update(kwargs)
        # We create initial string that we will gradually extend while we
        # process the data and in the end it will be written into a pdb file.
        string = "REMARK File generated using pyWINDOW."
        # Number of items (atoms) in the provided system.
        len_ = system[settings['atom_ids']].shape[0]
        # We process the remarks, if any, given by the user (optional).
        if isinstance(settings['remarks'], (list, tuple)):
            # If a list or tuple of remarks each is written at a new line
            # with the REMARK prefix not to have to long remark line.
            for remark in settings['remarks']:
                string = "\n".join([string, 'REMARK {0}'.format(remark)])
        else:
            # Otherwise if it's a single string or an int/float we just write
            # it under single remark line, otherwise nothing happens.
            if isinstance(settings['remarks'], (str, int, float)):
                remark = settings['remarks']
                string = "\n".join([string, 'REMARK {0}'.format(remark)])
        # If there is a unit cell (crystal data) provided we need to add it.
        if settings['cryst'] in system.keys():
            if system[settings['cryst']].any():
                cryst_line = "CRYST1"
                cryst = system[settings['cryst']]
                # The user have to provide the crystal data as a list/array
                # of six items containing unit cell edges lengths a, b and c
                # in x, y and z directions and three angles, or it can be.
                # Other options are not allowed for simplicity. It can convert
                # from the lattice array using function from utilities.
                for i in cryst[:3]:
                    cryst_line = "".join([cryst_line, "{0:9.3f}".format(i)])
                for i in cryst[3:]:
                    cryst_line = "".join([cryst_line, "{0:7.2f}".format(i)])
                # This is kind of messy, by default the data written in PDB
                # file should be P1 symmetry group therefore containing all
                # atom coordinates and not considering symmetry operations.
                # But, user can still define a space group if he wishes to.
                if settings['space_group'] is not None:
                    space_group = settings['space_group']
                else:
                    space_group = "{0}".format("P1")
                cryst_line = " ".join([cryst_line, space_group])
                # We add the unit cell parameters to the main string.
                string = "\n".join([string, cryst_line])
        # For the sake of code readability we extract interesting data from the
        # system. Atom_ids are the atom ids written at the third column of a
        # PDB file and the user has here the freedom to use the forcefield
        # assigned ones. However, they have to specify it directly using the
        # atom_ids key. Otherwise, the 'elements' array from system object
        # will be used, that is also used for elements in the last column of
        # a PDB file. Other parameters like residue name (resName), chain id
        # (chainID) and residue sequence (resSeq) can be controlled by
        # appropriate parameter keyword passed to this function, Otherwise
        # the default values from settings dictionary are used.
        atom_ids = system[settings['atom_ids']]
        elements = system[settings['elements']]
        # If the 'elements' array of the system need deciphering atom keys this
        # is done if the user sets decipher to True. They can also provided
        # forcefield, otherwise it's None which equals to DLF.
        if settings['decipher'] is True:
            elements = np.array([
                decipher_atom_key(
                    key, forcefield=settings['forcefield']) for key in elements
            ])
        coordinates = system[settings['coordinates']]
        for i in range(len_):
            atom_line = "ATOM  {0:5d}".format(i + 1)
            atom_id = "{0:4}".format(atom_ids[i].center(4))
            resName = "{0:3}".format(settings['resName'])
            chainID = settings['chainID']
            atom_line = " ".join([atom_line, atom_id, resName, chainID])
            resSeq = str(settings['resSeq']).rjust(4)
            atom_line = "".join([atom_line, resSeq])
            coor = "{0:8.3f}{1:8.3f}{2:8.3f}".format(
                coordinates[i][0],
                coordinates[i][1],
                coordinates[i][2], )
            atom_line = "    ".join([atom_line, coor])
            big_space = "{0}".format(" ".center(22))
            element = "{0:2}  ".format(elements[i].rjust(2))
            atom_line = "".join([atom_line, big_space, element])
            string = "\n".join([string, atom_line])
        # The connectivity part is to be written after a function calculating
        # connectivity is finished
        # "Everything that has a beginning has an end" by Neo. :)
        string = "\n".join([string, 'END'])
        # Check if .pdb extension is missing from filepath.
        if filepath[-4:].lower() != '.pdb':
            filepath = ".".join((filepath, 'pdb'))
        # Write the string to a a PDB file.
        with open(filepath, 'w+') as file:
            file.write(string)
