import numpy as np
from transition1x import Dataloader
from ase.io import read
import h5py
from pymatgen.core.structure import Molecule


def load_data_from_h5(path_to_h5_file: str, ref_coords: np.ndarray, ref_atomic_nums: np.ndarray, ref_molecule: Molecule) -> list:
    """
    Load data from an HDF5 file, including positions, atomic numbers, and reference coordinates.

    Parameters
    ----------
    path_to_h5_file : str
        The path to the HDF5 file containing molecular data.
    ref_coords : np.ndarray
        Reference coordinates for molecule comparison.
    ref_atomic_nums : np.ndarray
        Reference atomic numbers for molecule comparison.
    ref_molecule : Molecule
        Reference molecule object for comparison.

    Returns
    -------
    list
        A list containing tuples of positions, atomic numbers, reference coordinates, reference atomic numbers,
        and reference molecule for each entry in the HDF5 file.

    Examples
    --------
    >>> data_list = load_data_from_h5('/path/to/file.h5', ref_coords, ref_atomic_nums, ref_molecule)
    """
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
        return []  # Return an empty list for files with errors


def load_reference_molecule(path_to_ref_molecule: str = None) -> tuple:
    """
    Load reference molecule information from an XYZ file.

    Parameters
    ----------
    path_to_ref_molecule : str, optional
        The path to the XYZ file containing the reference molecule. If not provided, returns None.

    Returns
    -------
    tuple
        A tuple containing reference coordinates, reference atomic numbers, and a reference Molecule object.

    Examples
    --------
    >>> reference_coords, reference_atomic_nums, reference_molecule = load_reference_molecule('/path/to/molecule.xyz')
    >>> print(reference_coords)
    >>> print(len(reference_coords))
    >>> print(reference_atomic_nums)
    """
    ase_object = read(path_to_ref_molecule)
    reference_coords = ase_object.get_positions()
    reference_atomic_nums = ase_object.get_atomic_numbers()
    reference_molecule = Molecule(reference_atomic_nums, reference_coords)
    return reference_coords, reference_atomic_nums, reference_molecule


if __name__ == '__main__':
    path_to_ref_molecule = '../../tests/264_noise00.xyz'
    reference_coords, reference_atomic_nums, reference_molecule = load_reference_molecule(path_to_ref_molecule)

    print(reference_coords)
    print(len(reference_coords))
    print(reference_atomic_nums)
