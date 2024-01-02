from pymatgen.core.structure import Molecule
from mol_graph_funcs import compare_mols


def process_molecules(data_batch):
    results = []

    # Load the reference molecule
    reference_molecule = load_reference_molecule()

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
