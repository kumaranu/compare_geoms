import os
import os.path
import numpy as np
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash


def compare_mols(molecule_1, molecule_2) -> bool:
    """
    Compare two molecules based on their graph structure.

    Args:
        molecule_1 (pymatgen.core.structure.Molecule): First molecule.
        molecule_2 (pymatgen.core.structure.Molecule): Second molecule.

    Returns:
        bool: True if the molecules have the same graph structure, False otherwise.
    """
    molgraph_1 = create_molecule_graph(molecule_1)
    molgraph_2 = create_molecule_graph(molecule_2)

    graph_1 = molgraph_1.graph.to_undirected()
    graph_2 = molgraph_2.graph.to_undirected()

    add_specie_suffix(graph_1)
    add_specie_suffix(graph_2)

    graph_1_hash = get_graph_hash(graph_1)
    graph_2_hash = get_graph_hash(graph_2)

    return graph_1_hash == graph_2_hash


def create_molecule_graph(molecule):
    """
    Create a molecule graph using the OpenBabelNN strategy.

    Args:
        molecule (pymatgen.core.structure.Molecule): The molecule.

    Returns:
        pymatgen.analysis.graphs.MoleculeGraph: The molecule graph.
    """
    return MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())


def add_specie_suffix(graph):
    """
    Add a suffix to each node's 'specie' attribute in the graph.

    Args:
        graph (networkx.Graph): The graph.
    """
    for idx in graph.nodes():
        graph.nodes()[idx]["specie"] = graph.nodes()[idx]["specie"] + str(idx)


def get_graph_hash(graph):
    """
    Get the hash of the graph using the Weisfeiler-Lehman algorithm.

    Args:
        graph (networkx.Graph): The graph.

    Returns:
        str: The graph hash.
    """
    return weisfeiler_lehman_graph_hash(graph, node_attr='specie')
