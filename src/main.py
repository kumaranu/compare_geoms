from multiprocessing import Pool
from molecule_processing import process_molecules
from data_loading import load_data_from_h5, load_reference_molecule
import time


def compare_one_ref_mol(path_to_h5_file: str, path_to_ref_molecule: str) -> None:
    """
    Compare the reference molecule to a set of molecules stored in an HDF5 file using multiprocessing.

    Parameters
    ----------
    path_to_h5_file : str
        The path to the HDF5 file containing the molecular data.
    path_to_ref_molecule : str
        The path to the reference molecule file (XYZ format).

    Returns
    -------
    None
        The function prints the results of the molecule comparison.

    Notes
    -----
    The function loads the reference molecule and molecular data, then uses multiprocessing to parallelize
    the molecule processing. Adjust the chunksize parameter for optimal performance.

    Examples
    --------
    >>> compare_one_ref_mol('path/to/molecules.h5', 'path/to/reference_molecule.xyz')
    """
    ref_coords, ref_atomic_nums, ref_molecule = load_reference_molecule(path_to_ref_molecule)
    data_list = load_data_from_h5(path_to_h5_file, ref_coords, ref_atomic_nums, ref_molecule)

    t1 = time.time()
    with Pool(10) as pool:
        chunksize = 100  # Adjust the chunksize based on experimentation
        results = pool.map(process_molecules, [data_list[i:i+chunksize] for i in range(0, len(data_list), chunksize)])
    print(f'{time.time()-t1} seconds')

    # Flatten the list of results
    results = [result for sublist in results for result in sublist]

    for result in results:
        print(result)


if __name__ == "__main__":
    """
    Run the comparison of a reference molecule to a set of molecules stored in an HDF5 file.

    This script initializes the paths to the HDF5 file and the reference molecule, then calls
    the compare_one_ref_mol function to perform the molecule comparison.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Examples
    --------
    >>> python script_name.py
    """
    path_to_ref_molecule = '/home/kumaranu/Documents/analysis/molecules_fromscratch_noised_renamed_b00/264_noise00.xyz'
    path_to_h5_file = '/home/kumaranu/Documents/compare_geoms/tests/output_9953.h5'
    compare_one_ref_mol(path_to_h5_file, path_to_ref_molecule)
