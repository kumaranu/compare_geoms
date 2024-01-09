from multiprocessing import Pool

import networkx as nx

from compare_geoms.molecule_processing import process_molecules
from compare_geoms.data_loading import load_data_from_h5, load_reference_molecule
import os, time
from pysmiles import read_smiles

from pymatgen.core.structure import Molecule
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.graphs import MoleculeGraph
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from pandas import read_csv

from mol_graph_funcs import create_molecule_graph,\
    add_specie_suffix, get_graph_hash


def compare_one_ref_mol_from_smiles(smiles_string: str, path_to_ref_molecule: str) -> bool:
    molecule_smiles = read_smiles(smiles_string, explicit_hydrogen=True,
                                  zero_order_bonds=True, reinterpret_aromatic=True)

    graph_1 = molecule_smiles.to_undirected()

    #for idx in graph_1.nodes():
    #    if 'element' in graph_1.nodes()[idx]:
    #        graph_1.nodes()[idx]['specie'] = graph_1.nodes()[idx].pop('element') + str(idx)

    ref_molecule = Molecule.from_file(path_to_ref_molecule)
    molgraph_2 = create_molecule_graph(ref_molecule)
    graph_2 = molgraph_2.graph.to_undirected()

    # for idx in graph_2.nodes():
    #     print(graph_2.nodes()[idx]["specie"])

    # print(graph_2.nodes['specie'])
    #add_specie_suffix(graph_2)
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
    return nx.is_isomorphic(graph_1, graph_2)


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


import logging
logging.getLogger('pysmiles').setLevel(logging.CRITICAL)  # Anything higher than warning


def process_molecule1(args):
    path_to_ref_molecule, smiles_strings = args
    for smiles_string in smiles_strings:
        if compare_one_ref_mol_from_smiles(smiles_string, path_to_ref_molecule):
            print(smiles_string, path_to_ref_molecule)
            return 1
    return 0


if __name__ == "__main__":
    path_to_h5_file = '../../tests/output_9953.h5'
    path_to_csv_file = '/home/kumaranu/Downloads/b97d3.csv'
    smiles_strings = read_csv(path_to_csv_file)['rsmi']

    inputs = []
    smiles_df = read_csv(path_to_csv_file)
    for i in range(265):
        path_to_ref_molecule = (f'/home/kumaranu/Documents/analysis/'
                                f'molecules_fromscratch_noised_renamed_b00/{i:03}_noise00.xyz')
        partial_process_molecule = inputs.append((path_to_ref_molecule, smiles_df['rsmi']))

    # Use multiprocessing with 50 pools
    with Pool(20) as pool:
            results = pool.map(process_molecule1, inputs)
