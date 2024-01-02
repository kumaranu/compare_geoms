from multiprocessing import Pool
from molecule_processing import process_molecules
from data_loading import load_data, load_reference_molecule


def main():
    # Load data and reference molecule
    data_list, reference_molecule = load_data()

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
