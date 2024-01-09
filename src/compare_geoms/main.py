from multiprocessing import Pool
from compare_geoms.molecule_processing import process_molecules
from compare_geoms.data_loading import load_data_from_h5, load_reference_molecule
import os, time
from pysmiles import read_smiles

from pymatgen.core.structure import Molecule
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.graphs import MoleculeGraph
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash

from mol_graph_funcs import create_molecule_graph,\
    add_specie_suffix, get_graph_hash


def compare_one_ref_mol_from_smiles(smiles_string: str, path_to_ref_molecule: str) -> bool:
    molecule_smiles = read_smiles(smiles_string, explicit_hydrogen=True,
                                  zero_order_bonds=True, reinterpret_aromatic=True)

    graph_1 = molecule_smiles.to_undirected()

    for idx in graph_1.nodes():
        if 'element' in graph_1.nodes()[idx]:
            graph_1.nodes()[idx]['specie'] = graph_1.nodes()[idx].pop('element') + str(idx)

    ref_molecule = Molecule.from_file(path_to_ref_molecule)
    molgraph_2 = create_molecule_graph(ref_molecule)
    graph_2 = molgraph_2.graph.to_undirected()

    # for idx in graph_2.nodes():
    #     print(graph_2.nodes()[idx]["specie"])

    # print(graph_2.nodes['specie'])
    add_specie_suffix(graph_2)
    # print('\n\n\n\n\naaaaaaaa\n\n\n\n')
    # for idx in graph_2.nodes():
    #     print(graph_2.nodes()[idx]["specie"])
    # print(graph_2)
    # print(graph_2.edges)
    graph_2_hash = get_graph_hash(graph_2)

    ref_molecule_graph = MoleculeGraph.with_local_env_strategy(ref_molecule, OpenBabelNN())

    # ref_molecule_graph = ref_molecule_graph.
    # graph_1, graph_2 = map(lambda graph: graph.graph.to_undirected(), [molecule_graph, ref_molecule_graph])
    # print(graph_1.nodes)
    # mol.nodes(data='element'))

    # for idx in :
    #     print(idx)
    #     # graph.nodes()[idx]["specie"] += str(idx)

    # Calculate and compare graph hashes using Weisfeiler-Lehman algorithm
    # graph_1_hash = weisfeiler_lehman_graph_hash(graph_1, node_attr='specie')
    # graph_2_hash = weisfeiler_lehman_graph_hash(graph_2, node_attr='specie')

    # return graph_1_hash == graph_2_hash
    return 0


def compare_one_ref_mol(path_to_data: str, path_to_ref_molecule: str) -> None:
    """
    Compare the reference molecule to a set of molecules stored in HDF5 file(s) using multiprocessing.

    Parameters
    ----------
    path_to_data : str
        The path to the HDF5 file or directory containing the molecular data.
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
    >>> compare_one_ref_mol('path/to/molecule_directory', 'path/to/reference_molecule.xyz')
    """
    # Check if the given path is a directory
    if os.path.isdir(path_to_data):
        h5_files = [os.path.join(path_to_data, file) for file in os.listdir(path_to_data) if file.endswith('.h5')]
    else:
        h5_files = [path_to_data]

    ref_coords, ref_atomic_nums, ref_molecule = load_reference_molecule(path_to_ref_molecule)

    for h5_file in h5_files:
        data_list = load_data_from_h5(h5_file, ref_coords, ref_atomic_nums, ref_molecule)

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
    path_to_ref_molecule = '../../tests/264_noise00.xyz'
    path_to_h5_file = '../../tests/output_9953.h5'
    # compare_one_ref_mol(path_to_h5_file, path_to_ref_molecule)
    smiles_string = '[C:1]([c:2]1[n:3][o:4][n:5][n:6]1)([H:7])([H:8])[H:9]'
    compare_one_ref_mol_from_smiles(smiles_string, path_to_ref_molecule)
