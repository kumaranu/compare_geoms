from transition1x import Dataloader
from pymatgen.core.structure import Molecule
from ase.io import read
from aa import compare_mols
from multiprocessing import Pool


def process_molecules(data_batch):
    results = []

    for data in data_batch:
        coords1 = data.get('positions', [])
        atomic_nums1 = data.get('atomic_numbers', [])

        # Check if the number of atoms is the same
        if len(atomic_nums1) != len(reference_atomic_nums):
            continue

        # Check if the atomic numbers are the same
        if set(atomic_nums1) != set(reference_atomic_nums):
            continue

        molecule1 = Molecule(atomic_nums1, coords1)

        # Compare each molecule with the reference molecule
        if compare_mols(molecule1, reference_molecule):
            results.append("same")

    return results


def main():
    path_to_h5_file = '/home/kumaranu/Documents/transition1x.h5'
    dataloader = Dataloader(path_to_h5_file)

    ase_object = read('molecules_fromscratch_noised_renamed_b00/264_noise00.xyz')
    reference_coords = ase_object.get_positions()
    global reference_atomic_nums
    reference_atomic_nums = ase_object.get_atomic_numbers()

    # Create a reference molecule
    global reference_molecule
    reference_molecule = Molecule(reference_atomic_nums, reference_coords)

    # Prepare the data for parallel processing
    data_list = list(dataloader)

    # Use multiprocessing to parallelize the loop
    with Pool(10) as pool:
        chunksize = 100  # Adjust the chunksize based on experimentation
        results = pool.map(process_molecules, [data_list[i:i+chunksize] for i in range(0, len(data_list), chunksize)])

    # Flatten the list of results
    results = [result for sublist in results for result in sublist]

    for result in results:
        print(result)

if __name__ == "__main__":
    main()

