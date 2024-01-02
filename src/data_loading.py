from transition1x import Dataloader
from ase.io import read


def load_data():
    path_to_h5_file = '/home/kumaranu/Documents/transition1x.h5'
    dataloader = Dataloader(path_to_h5_file)

    ase_object = read('molecules_fromscratch_noised_renamed_b00/264_noise00.xyz')
    reference_coords = ase_object.get_positions()
    reference_atomic_nums = ase_object.get_atomic_numbers()

    # Create a reference molecule
    reference_molecule = Molecule(reference_atomic_nums, reference_coords)

    # Load data
    data_list = list(dataloader)

    return data_list, reference_molecule


def load_reference_molecule():
    ase_object = read('molecules_fromscratch_noised_renamed_b00/264_noise00.xyz')
    reference_coords = ase_object.get_positions()
    reference_atomic_nums = ase_object.get_atomic_numbers()

    # Create a reference molecule
    reference_molecule = Molecule(reference_atomic_nums, reference_coords)

    return reference_molecule
