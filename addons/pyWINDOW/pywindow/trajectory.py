#!/usr/bin/env python3
import numpy as np
from mmap import mmap, ACCESS_READ
from contextlib import closing
from multiprocessing import Pool

from .io_tools import Input, Output
from .utilities import is_number, create_supercell, lattice_matrix_to_unit_cell
from .molecular import MolecularSystem


class _ParallelAnalysisError(Exception):
    def __init__(self, message):
        self.message = message


def make_supercell(system, matrix, supercell=[1, 1, 1]):
    user_supercell = [[1, supercell[0]], [1, supercell[1]], [1, supercell[1]]]
    system = create_supercell(system, matrix, supercell=user_supercell)
    return MolecularSystem.load_system(system)


class DLPOLY(object):
    def __init__(self):
        self.frames = {}
        self.analysis_output = {}
        # Image conversion - periodic boundary key.
        self._imcon = {
            0: 'nonperiodic',
            1: 'cubic',
            2: 'orthorhombic',
            3: 'parallelepiped',
            4: 'truncated octahedral',
            5: 'rhombic dodecahedral',
            6: 'x-y parallelogram',
            7: 'hexagonal prism',
        }
        # Trajectory key - content type.
        self._keytrj = {
            0: 'coordinates',
            1: 'coordinates and velocities',
            2: 'coordinates, velocities and forces',
        }

    def load_HISTORY(self, file_path):
        """ Return file_path as a class attribute"""
        self.file_path = file_path
        self.trajectory_map = {}

        with open(file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                progress = 0
                line = 0
                frame = 0
                cell_param_line = 0
                # We need to first process trajectory file's header.
                header_flag = True
                while progress <= len(mapped_file):
                    line = line + 1
                    # We read a binary data from a mapped file.
                    bline = mapped_file.readline()
                    # If the bline length equals zero we terminate.
                    # We reached end of the file but still add the last frame!
                    if len(bline) == 0:
                        self.trajectory_map[frame] = [frame_start, progress]
                        frame = frame + 1
                        break
                    # We need to decode byte line into an utf-8 string.
                    sline = bline.decode("utf-8").strip('\n').split()
                    # We extract map's byte coordinates for each frame
                    if header_flag is False:
                        if sline[0] == 'timestep':
                            self.trajectory_map[frame] = [
                                frame_start, progress
                            ]
                            frame_start = progress
                            frame = frame + 1
                    # Here we extract the map's byte coordinates for the header
                    # And also the periodic system type needed for later.
                    if header_flag is True:
                        if sline[0] == 'timestep':
                            self.trajectory_map['header'] = self._decode_head(
                                [0, progress])
                            frame_start = progress
                            header_flag = False
                    progress = progress + len(bline)
            self.no_of_frames = frame
            self.get_frames = self._get_frames
        return self

    def _decode_head(self, header_coordinates):
        start, end = header_coordinates
        with open(self.file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                header = [
                    i.split()
                    for i in mapped_file[start:end].decode("utf-8").split('\n')
                ]
                header = [int(i) for i in header[1]]
        self.periodic_boundary = self._imcon[header[1]]
        self.content_type = self._keytrj[header[0]]
        self.no_of_atoms = header[2]
        return header

    def _get_frames(self, frames, extract_data=True, override=False):
        if override is True:
            self.frames = {}
        if isinstance(frames, int):
            frame = self._get_frame(self.trajectory_map[frames], extract_data)
            if frames not in self.frames.keys():
                self.frames[frames] = frame
            return frame
        if isinstance(frames, list):
            for frame in frames:
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, tuple):
            for frame in range(frames[0], frames[1]):
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, str):
            if frames in ['all', 'everything']:
                for frame in range(0, self.no_of_frames):
                    if frame not in self.frames.keys():
                        self.frames[frame] = self._get_frame(
                            self.trajectory_map[frame], extract_data)

    def _get_frame(self, frame_coordinates, extract_data):
        start, end = frame_coordinates
        with open(self.file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                if extract_data is False:
                    return mapped_file[start:end].decode("utf-8")
                else:
                    # [:-1] because the split results in last list empty.
                    frame = [
                        i.split()
                        for i in mapped_file[start:end].decode("utf-8").split(
                            '\n')
                    ][:-1]
                    return self._decode_frame(frame)

    def _decode_frame(self, frame):
        frame_data = {
            'frame_info': {
                'nstep': int(frame[0][1]),
                'natms': int(frame[0][2]),
                'keytrj': int(frame[0][3]),
                'imcon': int(frame[0][4]),
                'tstep': float(frame[0][5])
            }
        }
        start_line = 1
        if frame_data['frame_info']['imcon'] in [1, 2, 3]:
            frame_data['lattice'] = np.array(frame[1:4], dtype=float)
            frame_data['unit_cell'] = lattice_matrix_to_unit_cell(frame_data[
                'lattice'])
            start_line = 4
        # Depending on what the trajectory key is (see __init__) we need
        # to extract every second/ third/ fourth line for elements and coor.
        elements = []
        coordinates = []
        velocities = []
        forces = []
        for i in range(len(frame[start_line:])):
            i_ = i + start_line
            if frame_data['frame_info']['keytrj'] == 0:
                if i % 2 == 0:
                    elements.append(frame[i_][0])
                if i % 2 == 1:
                    coordinates.append(frame[i_])
            if frame_data['frame_info']['keytrj'] == 1:
                if i % 3 == 0:
                    elements.append(frame[i_][0])
                if i % 3 == 1:
                    coordinates.append(frame[i_])
                if i % 3 == 2:
                    velocities.append(frame[i_])
            if frame_data['frame_info']['keytrj'] == 2:
                if i % 4 == 0:
                    elements.append(frame[i_][0])
                if i % 4 == 1:
                    coordinates.append(frame[i_])
                if i % 4 == 2:
                    velocities.append(frame[i_])
                if i % 4 == 3:
                    forces.append(frame[i_])
        frame_data['atom_ids'] = np.array(elements)
        frame_data['coordinates'] = np.array(coordinates, dtype=float)
        if velocities:
            frame_data['velocities'] = np.array(velocities, dtype=float)
        if forces:
            frame_data['forces'] = np.array(forces, dtype=float)
        return MolecularSystem.load_system(frame_data)

    def analysis(self, frames='all', ncpus=1, override=False, **kwargs):
        if override is True:
            self.analysis_output = {}
        if isinstance(frames, int):
            analysed_frame = self._analysis_serial(frames, ncpus, **kwargs)
            if frames not in self.analysis_output.keys():
                self.analysis_output[frames] = analysed_frame
            return analysed_frame
        else:
            frames_for_analysis = []
            if isinstance(frames, list):
                for frame in frames:
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, tuple):
                for frame in range(frames[0], frames[1]):
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, str):
                if frames in ['all', 'everything']:
                    for frame in range(0, self.no_of_frames):
                        if frame not in self.analysis_output.keys():
                            frames_for_analysis.append(frame)
            self._analysis_parallel(frames_for_analysis, ncpus, **kwargs)

    def _analysis_serial(self, frame, ncpus, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(
                            ncpus=ncpus, **kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
            else:
                results[molecule] = mol.full_analysis(ncpus=ncpus, **kwargs)
        return results

    def _analysis_parallel_execute(self, frame, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(**kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(**kwargs)
            else:
                results[molecule] = mol.full_analysis(**kwargs)
        return frame, results

    def _analysis_parallel(self, frames, ncpus, **kwargs):
        try:
            pool = Pool(processes=ncpus)
            parallel = [
                pool.apply_async(
                    self._analysis_parallel_execute,
                    args=(frame, ),
                    kwds=kwargs) for frame in frames
            ]
            results = [p.get() for p in parallel if p.get()[1] is not None]
            pool.terminate()
            for i in results:
                self.analysis_output[i[0]] = i[1]
        except TypeError:
            pool.terminate()
            raise _ParallelAnalysisError("Parallel analysis failed.")

    def check_HISTORY_for_errors(self, file_path):
        """
        """

        line = 0
        binary_step = 0
        timestep = 0
        timestep_flag = 'timestep'
        frame = 0
        warnings = []
        errors = []
        progress = 0

        warning_1 = "Warning 1: No comment line is present as the file header."
        warning_2 = "Warning 2: Second header line is missing from the file "
        warning_2 += "that contains information on the system's periodicity "
        warning_2 += "and the type of the trajectory file."
        warning_3 = "Warning 3: Comment line encountered in the middle of "
        warning_3 += "the trajectory file."

        error_1 = "Error 1: The trajectory is discontinous."
        error_2 = "Error 2: The file contains an empty line."

        with open(self.file_path, 'r') as trajectory_file:
            # We open the HISTORY trajectory file
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as file_binary_map:
                # We use this binary mapping feature that instead of loading
                # the full file into memory beforehand it only
                # maps the content. Especially useful with enormous files
                while binary_step < len(file_binary_map):
                    line += 1
                    binary_line = file_binary_map.readline()
                    binary_step = binary_step + len(binary_line)
                    progress_old = progress
                    progress = round(binary_step * 100 / len(file_binary_map),
                                     0)
                    string_line = binary_line.decode("utf-8").strip(
                        '\n').split()

                    # Warning 1
                    if line == 1:
                        if string_line[0] != 'DLFIELD':
                            warnings.append("Line {0}. {1}".format(line,
                                                                   warning_1))

                    # Warning 2
                    if line == 2:
                        if len(string_line) != 3:
                            warnings.append("Line {0}. {1}".format(line,
                                                                   warning_2))

                    # Error 1
                    if string_line:
                        if string_line[0] == timestep_flag:
                            frame += 1
                            old_timestep = timestep
                            timestep = int(string_line[1])
                            if old_timestep > timestep:
                                errors.append("Line {0}. {1}".format(line,
                                                                     error_1))

                    # Error 2
                    if len(string_line) == 0:
                        errors.append("Line {0}. {1}".format(line, error_2))

        if len(warnings) == 0 and len(errors) == 0:
            print('\nTrajectory file check finished. No warnings.')
            print('Frames: {0}'.format(frame))

        else:
            print('\nTrajectory file check finished.')
            print('Frames: {0}'.format(frame))
            for i in warnings:
                print(i)
            for i in errors:
                print(i)
            user_guide = "For more detailed description of warnings and "
            user_guide += "errors please see the pyWindow User's Guide."
            print(user_guide)
        return self


class XYZ(object):
    def __init__(self):
        self.frames = {}
        self.analysis_output = {}

    def load(self, file_path):
        """ Return file_path as a class attribute"""
        self.file_path = file_path
        self.trajectory_map = {}
        with open(file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                progress = 0
                line = 0
                frame = -1
                frame_start = 0
                while progress <= len(mapped_file):
                    line = line + 1
                    # We read a binary data from a mapped file.
                    bline = mapped_file.readline()
                    # If the bline length equals zero we terminate.
                    # We reached end of the file but still add the last frame!
                    if len(bline) == 0:
                        frame = frame + 1
                        self.trajectory_map[frame] = [frame_start, progress]
                        break
                    # We need to decode byte line into an utf-8 string.
                    sline = bline.decode("utf-8").strip('\n').split()
                    # We extract map's byte coordinates for each frame
                    if (len(sline) == 1 and is_number(sline[0]) and
                            progress > 0):
                        frame = frame + 1
                        self.trajectory_map[frame] = [frame_start, progress]
                        frame_start = progress
                    # Here we extract the map's byte coordinates for the header
                    # And also the periodic system type needed for later.
                    progress = progress + len(bline)
            self.no_of_frames = frame + 1
            self.get_frames = self._get_frames
        return self

    def _get_frames(self, frames, extract_data=True, override=False):
        if override is True:
            self.frames = {}
        if isinstance(frames, int):
            frame = self._get_frame(self.trajectory_map[frames], extract_data)
            if frames not in self.frames.keys():
                self.frames[frames] = frame
            return frame
        if isinstance(frames, list):
            for frame in frames:
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, tuple):
            for frame in range(frames[0], frames[1]):
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, str):
            if frames in ['all', 'everything']:
                for frame in range(0, self.no_of_frames):
                    if frame not in self.frames.keys():
                        self.frames[frame] = self._get_frame(
                            self.trajectory_map[frame], extract_data)

    def _get_frame(self, frame_coordinates, extract_data):
        start, end = frame_coordinates
        with open(self.file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                if extract_data is False:
                    return mapped_file[start:end].decode("utf-8")
                else:
                    # [:-1] because the split results in last list empty.
                    frame = [
                        i.split()
                        for i in mapped_file[start:end].decode("utf-8").split(
                            '\n')
                    ][:-1]
                    return self._decode_frame(frame)

    def _decode_frame(self, frame):
        frame_data = {
            'frame_info': {
                'natms': int(frame[0][0]),
                'remarks': " ".join([*frame[1]]),
            }
        }
        start_line = 2
        elements = []
        coordinates = []
        for i in range(start_line, len(frame)):
            elements.append(frame[i][0])
            coordinates.append(frame[i][1:])
        frame_data['atom_ids'] = np.array(elements)
        frame_data['coordinates'] = np.array(coordinates, dtype=float)
        return MolecularSystem.load_system(frame_data)

    def analysis(self, frames='all', ncpus=1, override=False, **kwargs):
        if override is True:
            self.analysis_output = {}
        if isinstance(frames, int):
            analysed_frame = self._analysis_serial(frames, ncpus, **kwargs)
            if frames not in self.analysis_output.keys():
                self.analysis_output[frames] = analysed_frame
            return analysed_frame
        else:
            frames_for_analysis = []
            if isinstance(frames, list):
                for frame in frames:
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, tuple):
                for frame in range(frames[0], frames[1]):
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, str):
                if frames in ['all', 'everything']:
                    for frame in range(0, self.no_of_frames):
                        if frame not in self.analysis_output.keys():
                            frames_for_analysis.append(frame)
            self._analysis_parallel(frames_for_analysis, ncpus, **kwargs)

    def _analysis_serial(self, frame, ncpus, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(
                            ncpus=ncpus, **kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
            else:
                results[molecule] = mol.full_analysis(ncpus=ncpus, **kwargs)
        return results

    def _analysis_parallel_execute(self, frame, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(**kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(**kwargs)
            else:
                results[molecule] = mol.full_analysis(**kwargs)
        return frame, results

    def _analysis_parallel(self, frames, ncpus, **kwargs):
        try:
            pool = Pool(processes=ncpus)
            parallel = [
                pool.apply_async(
                    self._analysis_parallel_execute,
                    args=(frame, ),
                    kwds=kwargs) for frame in frames
            ]
            results = [p.get() for p in parallel if p.get()[1] is not None]
            pool.terminate()
            for i in results:
                self.analysis_output[i[0]] = i[1]
        except TypeError:
            pool.terminate()
            raise _ParallelAnalysisError("Parallel analysis failed.")


class PDB(object):
    def __init__(self):
        self.frames = {}
        self.analysis_output = {}

    def load(self, file_path):
        """ Return file_path as a class attribute"""
        self.file_path = file_path
        self.trajectory_map = {}
        with open(file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                progress = 0
                line = 0
                frame = -1
                frame_start = 0
                while progress <= len(mapped_file):
                    line = line + 1
                    # We read a binary data from a mapped file.
                    bline = mapped_file.readline()
                    # If the bline length equals zero we terminate.
                    # We reached end of the file but still add the last frame!
                    if len(bline) == 0:
                        frame = frame + 1
                        if progress - frame_start > 10:
                            self.trajectory_map[frame] = [
                                frame_start, progress
                            ]
                        break
                    # We need to decode byte line into an utf-8 string.
                    sline = bline.decode("utf-8").strip('\n').split()
                    # We extract map's byte coordinates for each frame
                    if len(sline) == 1 and sline[0] == 'END':
                        frame = frame + 1
                        self.trajectory_map[frame] = [frame_start, progress]
                        frame_start = progress
                    # Here we extract the map's byte coordinates for the header
                    # And also the periodic system type needed for later.
                    progress = progress + len(bline)
            self.no_of_frames = frame
            self.get_frames = self._get_frames
        return self

    def _get_frames(self, frames, extract_data=True, override=False):
        if override is True:
            self.frames = {}
        if isinstance(frames, int):
            frame = self._get_frame(self.trajectory_map[frames], extract_data)
            if frames not in self.frames.keys():
                self.frames[frames] = frame
            return frame
        if isinstance(frames, list):
            for frame in frames:
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, tuple):
            for frame in range(frames[0], frames[1]):
                if frame not in self.frames.keys():
                    self.frames[frame] = self._get_frame(
                        self.trajectory_map[frame], extract_data)
        if isinstance(frames, str):
            if frames in ['all', 'everything']:
                for frame in range(0, self.no_of_frames):
                    if frame not in self.frames.keys():
                        self.frames[frame] = self._get_frame(
                            self.trajectory_map[frame], extract_data)

    def _get_frame(self, frame_coordinates, extract_data):
        start, end = frame_coordinates
        with open(self.file_path, 'r') as trajectory_file:
            with closing(
                    mmap(
                        trajectory_file.fileno(), 0,
                        access=ACCESS_READ)) as mapped_file:
                if extract_data is False:
                    return mapped_file[start:end].decode("utf-8")
                else:
                    # In case of PDB we do not split lines!
                    frame = mapped_file[start:end].decode("utf-8").split('\n')
                    return self._decode_frame(frame)

    def _decode_frame(self, frame):
        frame_data = {}
        elements = []
        coordinates = []
        for i in range(len(frame)):
            if frame[i][:6] == 'REMARK':
                if 'REMARKS' not in frame_data.keys():
                    frame_data['REMARKS'] = []
                frame_data['REMARKS'].append(frame[i][6:])
            if frame[i][:6] == 'CRYST1':
                cryst = np.array(
                    [
                        frame[i][6:15], frame[i][15:24], frame[i][24:33],
                        frame[i][33:40], frame[i][40:47], frame[i][47:54]
                    ],
                    dtype=float)
                # This is in case of nonperiodic systems, often they have
                # a,b,c unit cell parameters as 0,0,0.
                if sum(cryst[0:3]) != 0:
                    frame_data['CRYST1'] = cryst
            if frame[i][:6] in ['HETATM', 'ATOM  ']:
                elements.append(frame[i][12:16].strip())
                coordinates.append(
                    [frame[i][30:38], frame[i][38:46], frame[i][46:54]])
        frame_data['atoms_ids'] = np.array(elements, dtype='<U8')
        frame_data['coordinates'] = np.array(coordinates, dtype=float)
        return MolecularSystem.load_system(frame_data)

    def analysis(self, frames='all', ncpus=1, override=False, **kwargs):
        if override is True:
            self.analysis_output = {}
        if isinstance(frames, int):
            analysed_frame = self._analysis_serial(frames, ncpus, **kwargs)
            if frames not in self.analysis_output.keys():
                self.analysis_output[frames] = analysed_frame
            return analysed_frame
        else:
            frames_for_analysis = []
            if isinstance(frames, list):
                for frame in frames:
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, tuple):
                for frame in range(frames[0], frames[1]):
                    if frame not in self.analysis_output.keys():
                        frames_for_analysis.append(frame)
            if isinstance(frames, str):
                if frames in ['all', 'everything']:
                    for frame in range(0, self.no_of_frames):
                        if frame not in self.analysis_output.keys():
                            frames_for_analysis.append(frame)
            self._analysis_parallel(frames_for_analysis, ncpus, **kwargs)

    def _analysis_serial(self, frame, ncpus, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(
                            ncpus=ncpus, **kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(
                                ncpus=ncpus, **kwargs)
            else:
                results[molecule] = mol.full_analysis(ncpus=ncpus, **kwargs)
        return results

    def _analysis_parallel_execute(self, frame, **kwargs):
        frame_ = self._get_frame(self.trajectory_map[frame], extract_data=True)
        molecular_system = MolecularSystem()
        molecular_system.load_trajectory_frame(frame_)
        if 'swap_atoms' in kwargs:
            molecular_system.swap_atom_keys(kwargs['swap_atoms'])
        if 'forcefield' in kwargs:
            molecular_system.decipher_atom_keys(kwargs['forcefield'])
        molecular_system.separate_discrete_molecules()
        results = {}
        for molecule in molecular_system.molecules:
            mol = molecular_system.molecules[molecule]
            if 'size' in kwargs:
                size = kwargs['size']
                if isinstance(size, int):
                    if mol.no_of_atoms == size:
                        results[molecule] = mol.full_analysis(**kwargs)
                if isinstance(size, tuple) and isinstance(size[0], str):
                    if size[0] in ['bigger', 'greater', 'larger', 'more']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['smaller', 'less']:
                        if mol.no_of_atoms > size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['not', 'isnot', 'notequal', 'different']:
                        if mol.no_of_atoms != size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['is', 'equal', 'exactly']:
                        if mol.no_of_atoms == size[1]:
                            results[molecule] = mol.full_analysis(**kwargs)
                    if size[0] in ['between', 'inbetween']:
                        if size[1] < mol.no_of_atoms < size[2]:
                            results[molecule] = mol.full_analysis(**kwargs)
            else:
                results[molecule] = mol.full_analysis(**kwargs)
        return frame, results

    def _analysis_parallel(self, frames, ncpus, **kwargs):
        try:
            pool = Pool(processes=ncpus)
            parallel = [
                pool.apply_async(
                    self._analysis_parallel_execute,
                    args=(frame, ),
                    kwds=kwargs) for frame in frames
            ]
            results = [p.get() for p in parallel if p.get()[1] is not None]
            pool.terminate()
            for i in results:
                self.analysis_output[i[0]] = i[1]
        except TypeError:
            pool.terminate()
            raise _ParallelAnalysisError("Parallel analysis failed.")
