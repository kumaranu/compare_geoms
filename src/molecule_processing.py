import numpy as np
from pymatgen.core.structure import Molecule
from mol_graph_funcs import compare_mols
from data_loading import load_data_from_h5, load_reference_molecule
from multiprocessing import Pool


def process_molecules(data_batch):
    results = []
    for data_tuple in data_batch:
        coords, atomic_nums, ref_coords, ref_atomic_nums, ref_molecule = data_tuple

        # Check if the number of atoms is the same
        if len(atomic_nums) != len(ref_atomic_nums):
            return results

        # Check if the atomic numbers are the same
        if set(atomic_nums) != set(ref_atomic_nums):
            return results

        molecule_from_dataset = Molecule(atomic_nums, coords)

        # Compare each molecule with the reference molecule
        if compare_mols(molecule_from_dataset, ref_molecule):
            results.append("same")

    return results


if __name__ == '__main__':
    path_to_ref_molecule = '/home/kumaranu/Documents/analysis/molecules_fromscratch_noised_renamed_b00/264_noise00.xyz'
    path_to_h5_file = '/home/kumaranu/Documents/compare_geoms/tests/output_9953.h5'
    ref_coords, ref_atomic_nums, ref_molecule = load_reference_molecule(path_to_ref_molecule)
    data_list = load_data_from_h5(path_to_h5_file, ref_coords, ref_atomic_nums, ref_molecule)

    with Pool(10) as pool:
        # results = pool.map(process_molecules, data_list)
        chunksize = 100  # Adjust the chunksize based on experimentation
        results = pool.map(process_molecules, [data_list[i:i+chunksize] for i in range(0, len(data_list), chunksize)])
