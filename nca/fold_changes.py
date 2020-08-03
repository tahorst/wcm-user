#! /usr/bin/env python

"""
Use Network Component Analysis (NCA) to determine TF regulatory impacts on
each gene based on expression in all conditions.
"""

import argparse
import csv
import os

import numpy as np
from typing import Dict, List

import nca


BASE_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')

# Output related
OUTPUT_DIR = os.path.join(BASE_DIR, 'out')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)
RESULTS_FILE_TEMPLATE = os.path.join(OUTPUT_DIR, '{}nca_results.tsv')

# RegulonDB related
REGULON_DB_DIR = os.path.join(DATA_DIR, 'regulon-db')
REGULON_DB_SCRIPT = os.path.join(BASE_DIR, 'download_regulondb.sh')
TF_GENE_FILE = 'tf_genes.tsv'

# Sequencing related
COMPENDIUM_DIR = os.path.join(DATA_DIR, 'compendium')
GENE_NAMES_FILE = os.path.join(COMPENDIUM_DIR, 'gene_names.tsv')
SEQ_DIR = os.path.join(COMPENDIUM_DIR, 'seq')
SEQ_FILES = [
    'EcoMAC.tsv',
    # 'RNAseqEcoMACFormat.tsv',
    # 'GSE29076.tsv',
    # 'GSE72525.tsv',
    # 'GSE55662.tsv',
    # 'GSE50529.tsv',
    # 'GSE55365aerobic.tsv',
    # 'GSE55365anaerobic.tsv',
    ]


def load_regulon_db_file(filename: str) -> List[List[str]]:
    """
    Perform check to see if regulonDB data needs to be downloaded before loading
    a tsv file from the regulonDB data directory.

    Args:
        filename: name of file in regulonDB data directory to load

    Returns:
        data: data from a loaded tsv file
    """

    # Check if data exists or if the script needs to be run
    path = os.path.join(REGULON_DB_DIR, filename)
    if not os.path.exists(path):
        raise IOError(f'{path} does not exist. Run {REGULON_DB_SCRIPT} to download regulonDB data.')

    # Load data
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        data = [line for line in reader if not line[0].startswith('#')]

    return data

def load_gene_names() -> np.ndarray:
    """
    Loads genes names associated with sequencing data rows.

    Returns:
        genes: names of each gene, ordered to match sequencing data rows
    """

    with open(GENE_NAMES_FILE) as f:
        reader = csv.reader(f, delimiter='\t')
        genes = np.array([line[2] for line in reader])

    return genes

def load_seq_data() -> np.ndarray:
    """
    Load sequencing data from the compendium.

    Returns:
        matrix of normalized sequencing data representing log2 counts (n genes, m samples)

    TODO:
    - normalize data to include other files instead of just EcoMAC in SEQ_FILES
    """

    seq_data = []
    for filename in SEQ_FILES:
        path = os.path.join(SEQ_DIR, filename)

        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            seq_data.append(list(reader))

    return np.hstack(seq_data).astype(np.float64)

def load_tf_gene_interactions(verbose: bool = True) -> Dict[str, Dict[str, str]]:
    """
    Load regulonDB TF-gene interactions.

    Args:
        verbose: If True, prints warnings about loaded data

    Returns:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
    """

    data = load_regulon_db_file(TF_GENE_FILE)

    tf_genes = {}
    tf_idx = 0
    gene_idx = 1
    dir_idx = 2
    evidence_idx = 4  # TODO: filter on evidence?
    for line in data:
        # Extract columns of interest
        tf = line[tf_idx]
        gene = line[gene_idx]
        effect = line[dir_idx]

        # Check type of regulation
        if effect == 'activator':
            direction = 1
        elif effect == 'repressor':
            direction = -1
        elif effect == 'unknown':
            direction = 1  # TODO: handle this differently?
        else:
            raise ValueError(f'Unknown TF effect: {effect}')

        # Store new data
        genes = tf_genes.get(tf, {})
        if verbose and genes.get(gene, direction) != direction:
            print(f'Inconsistent regulatory effects for {tf} on {gene}')  # TODO: handle?
        genes[gene] = direction
        tf_genes[tf] = genes

    return tf_genes

def create_tf_map(
        gene_names: np.ndarray,
        tf_genes: Dict[str, Dict[str, str]],
        verbose: bool = True,
        ) -> (np.ndarray, np.ndarray):
    """
    Create an initial map reflecting the known regulatory network topology between
    transcription factors and genes. Also includes a column of ones for constitutive
    expression of each gene.

    Args:
        gene_names: gene IDs corresponding to each row of expression matrix
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        verbose: If True, prints warnings about loaded data

    Returns:
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
        tfs: IDs of TFs associated with each column of mapping
    """

    np.random.seed(0)

    gene_idx = {name: idx for idx, name in enumerate(gene_names)}
    n_genes = len(gene_names)
    n_tfs = len(tf_genes)
    mapping = np.zeros((n_genes, n_tfs))

    # Populate matrix with links between TFs and genes
    tfs = []
    for j, (tf, genes) in enumerate(sorted(tf_genes.items())):
        for gene, direction in genes.items():
            if gene not in gene_idx:
                if verbose:
                    print(f'Unknown gene: {gene}')
                continue

            mapping[gene_idx[gene], j] = np.random.rand() + 1
            # mapping[gene_idx[gene], j] = direction  # TODO: does this need a sign
        tfs.append(tf)

    return mapping, np.array(tfs)

def save_regulation(
        A: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        filename: str,
        ) -> None:
    """
    Save the results of NCA to a file showing each TF/gene pair.

    Args:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        filename: path to output file to save data to
    """

    relation_idx = np.where(A.T)
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['TF', 'Gene', 'FC'])
        for tf_idx, gene_idx in zip(*relation_idx):
            writer.writerow([tfs[tf_idx], genes[gene_idx], A[gene_idx, tf_idx]])

def parse_args() -> argparse.Namespace:
    """Parse command line args for options to run."""

    parser = argparse.ArgumentParser()

    default_nca = nca.METHODS[0]

    parser.add_argument('-f', '--force',
        action='store_true',
        help='If set, force a rerun of identifying the initial TF map, otherwise use cached values.')
    parser.add_argument('-l', '--label',
        help='Label prepended to output files for identification.')
    parser.add_argument('-m', '--method',
        choices=nca.METHODS,
        default=default_nca,
        help=f'NCA method to use, defined in nca.py (default: {default_nca}).')
    parser.add_argument('-v', '--verbose',
        action='store_true',
        help='If set, prints status updates for creating initial network map.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    # TODO: normalize seq data or only use EcoMAC data for now
    print('Loading data from files...')
    seq_data = load_seq_data()
    genes = load_gene_names()
    tf_genes = load_tf_gene_interactions(verbose=args.verbose)

    no_cache = True  # TODO: check if files exist
    if args.force or no_cache:
        print('Creating initial network mapping...')
        initial_tf_map, tfs = create_tf_map(genes, tf_genes, verbose=args.verbose)
        initial_tf_map, tfs = nca.nca_criteria_check(initial_tf_map, tfs, verbose=args.verbose)
        # TODO: save cache
    else:
        print('Loading cached initial network mapping...')
        # TODO: load cache
        raise NotImplementedError('Need to save and load cache.')

    A, P = getattr(nca, args.method)(seq_data, initial_tf_map)
    output_file = RESULTS_FILE_TEMPLATE.format(f'{args.label}_' if args.label else '')
    save_regulation(A, genes, tfs, output_file)
