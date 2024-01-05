import numpy as np
from transition1x import Dataloader
from ase.io import read
import h5py
from pymatgen.core.structure import Molecule

def load_data_from_h5(path_to_h5_file, ref_coords, ref_atomic_nums, ref_molecule):
    data_list = []
    try:
        with h5py.File(path_to_h5_file, 'r') as f:
            for i in f:
                positions = np.asarray(f[i]['positions'])
                atomic_numbers = np.asarray(f[i]['atomic_numbers'])
                data_list.append((positions, atomic_numbers, ref_coords, ref_atomic_nums, ref_molecule))
        return data_list
    except Exception as e:
        print(f"Error processing file {path_to_h5_file}: {e}")
        return []  # Return 0 for files with errors


def load_reference_molecule(path_to_ref_molecule=None):
    ase_object = read(path_to_ref_molecule)
    reference_coords = ase_object.get_positions()
    reference_atomic_nums = ase_object.get_atomic_numbers()
    reference_molecule = Molecule(reference_atomic_nums, reference_coords)
    return reference_coords, reference_atomic_nums, reference_molecule


if __name__ == '__main__':
    path_to_ref_molecule = '/home/kumaranu/Documents/analysis/molecules_fromscratch_noised_renamed_b00/264_noise00.xyz'
    reference_coords, reference_atomic_nums, reference_molecule = load_reference_molecule(path_to_ref_molecule)

    print(reference_coords)
    print(len(reference_coords))
    print(reference_atomic_nums)
