#! /usr/bin/env python

"""
Use Network Component Analysis (NCA) to determine TF regulatory impacts on
each gene based on expression in all conditions.
"""

import argparse
import csv
import os
import time
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

import nca


BASE_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')

# Output related
OUTPUT_DIR = os.path.join(BASE_DIR, 'out')
if not os.path.exists(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)
REGULATION_FILE = 'regulation.tsv'
ACTIVITY_FILE = 'activity.tsv'
HISTOGRAM_FILE = 'histogram.png'

# Cached results
NETWORK_CACHE_FILE = 'network.npy'
TF_CACHE_FILE = 'tfs.npy'

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

def load_seq_data(linearize: bool) -> np.ndarray:
    """
    Load sequencing data from the compendium.

    Args:
        linearize: if set, counts will be transformed from log2 space to linear space

    Returns:
        data: matrix of normalized sequencing data representing counts (n genes, m samples)
            linear if linearize, log2 otherwise

    TODO:
    - normalize data to include other files instead of just EcoMAC in SEQ_FILES
    """

    seq_data = []
    for filename in SEQ_FILES:
        path = os.path.join(SEQ_DIR, filename)

        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            seq_data.append(list(reader))

    data = np.hstack(seq_data).astype(np.float64)
    if linearize:
        data = 2**data

    return data

def load_tf_gene_interactions(
        split: bool = False,
        verbose: bool = True,
        ) -> Dict[str, Dict[str, int]]:
    """
    Load regulonDB TF-gene interactions.

    Args:
        split: if True, splits TFs into activator and repressor forms
        verbose: if True, prints warnings about loaded data

    Returns:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
            regulatory direction is 1 for positive regulation
            regulatory direction is -1 for negative regulation
            regulatory direction is 0 for unknown or conflicting regulation
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
            direction = 0
        else:
            raise ValueError(f'Unrecognized TF effect: {effect}')

        # Store new data
        if split:
            if effect == 'activator' or effect == 'unknown':
                split_tf = f'{tf}-activator'
                genes = tf_genes.get(split_tf, {})
                genes[gene] = 1
                tf_genes[split_tf] = genes

            if effect == 'repressor' or effect == 'unknown':
                split_tf = f'{tf}-repressor'
                genes = tf_genes.get(split_tf, {})
                genes[gene] = -1
                tf_genes[split_tf] = genes
        else:
            genes = tf_genes.get(tf, {})
            if genes.get(gene, direction) != direction:
                if verbose:
                    print(f'Inconsistent regulatory effects for {tf} on {gene}')
                genes[gene] = 0
            else:
                genes[gene] = direction
            tf_genes[tf] = genes

    return tf_genes

def create_tf_map(
        gene_names: np.ndarray,
        tf_genes: Dict[str, Dict[str, int]],
        verbose: bool = True,
        ) -> (np.ndarray, np.ndarray):
    """
    Create an initial map reflecting the known regulatory network topology between
    transcription factors and genes. Also includes a column of ones for constitutive
    expression of each gene.

    Args:
        gene_names: gene IDs corresponding to each row of expression matrix
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        verbose: if True, prints warnings about loaded data

    Returns:
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            positive regulation: positive number
            negative regulation: negative number
            unknown/ambiguous regulation: NaN
            no regulation: 0
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

            if direction == 0:
                mapping[gene_idx[gene], j] = np.nan
            else:
                mapping[gene_idx[gene], j] = direction * np.random.rand()
        tfs.append(tf)

    return mapping, np.array(tfs)

def save_regulation(
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        output_dir: str,
        ) -> None:
    """
    Save the results of NCA to a file showing each TF/gene pair.

    Args:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: TF activity for each condition
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        output_dir: path to directory to save the output
    """

    # Regulation data (TF-gene pairs)
    relation_idx = np.where(A.T)
    with open(os.path.join(output_dir, REGULATION_FILE), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['TF', 'Gene', 'FC'])
        for tf_idx, gene_idx in zip(*relation_idx):
            writer.writerow([tfs[tf_idx], genes[gene_idx], A[gene_idx, tf_idx]])

    # Activity data (TF in each condition)
    with open(os.path.join(output_dir, ACTIVITY_FILE), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Condition:'] + list(range(P.shape[1])))
        for tf, activity in zip(tfs, P):
            writer.writerow([tf] + list(activity))

def load_regulation(directory: str, genes: np.ndarray, tfs: np.ndarray) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Load the results of a previous NCA run saved to file with save_regulation.
    Return values should match the arguments passed to save_regulation.

    Args:
        directory: path to folder containing saved NCA results
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)

    Returns:
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: TF activity for each condition
    """

    with open(os.path.join(directory, REGULATION_FILE)) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        data = np.array(list(reader))

    fc_genes = data[:, headers.index('Gene')]
    fc_tfs = data[:, headers.index('TF')]
    fcs = data[:, headers.index('FC')].astype(float)

    # Recreate A matrix
    A = np.zeros((len(genes), len(tfs)))
    gene_idx = {gene: idx for idx, gene in enumerate(genes)}
    tf_idx = {tf: idx for idx, tf in enumerate(tfs)}
    for gene, tf, fc in zip(fc_genes, fc_tfs, fcs):
        A[gene_idx[gene], tf_idx[tf]] = fc

    # Recreate P matrix
    with open(os.path.join(directory, ACTIVITY_FILE)) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # strip headers
        P = np.array(list(reader))[:, 1:].astype(float)

    return A, P

def add_global_expression(
        tfs: np.ndarray,
        mapping: Optional[np.ndarray] = None
        ) -> (np.ndarray, np.ndarray):
    """
    Expand out TFs to include global expression (consituitive and regulated).
    This will capture genes not captured by TFs.

    Args:
        tfs: IDs of TFs associated with each column of mapping
        mapping: matrix representing network links between TFs and genes (n genes, m TFs)
            if None, no update is performed

    Returns:
        tfs: updated IDs with global expression IDs
        mapping: updated matrix with global expression columns

    TODO:
        - add sigma factors
        - add ppGpp
        - add noisy elements
    """

    tfs = np.hstack((np.array(['constituitive', 'regulated']), tfs))

    if mapping is not None:
        n_genes = mapping.shape[0]
        no_regulation = np.sum(mapping, axis=1) == 0
        constituitive = np.zeros((n_genes, 1))
        constituitive[no_regulation] = 1
        regulated = np.zeros((n_genes, 1))
        regulated[~no_regulation] = 1
        mapping = np.hstack((constituitive, regulated, mapping))

    return tfs, mapping

def match_statistics(
        E: np.ndarray,
        A: np.ndarray,
        P: np.ndarray,
        tf_genes: Dict[str, Dict[str, int]],
        genes: np.ndarray,
        tfs: np.ndarray,
        ) -> None:
    """
    Assess the accuracy of NCA results by printing statistics about how well
    the regulation direction and overall expression is captured.

    Args:
        E: original expression data (n genes, m conditions)
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF activity in each condition (o TFs, m conditions)
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)

    TODO:
        - stats for each TF
        - track skipped TFs and genes needed to satisfy NCA constraints
    """

    # Variable to track stats
    total_neg = 0
    correct_neg = 0
    dropped_neg = 0
    total_pos = 0
    correct_pos = 0
    dropped_pos = 0
    ambiguous = 0

    gene_idx = {gene: i for i, gene in enumerate(genes)}
    tf_idx = {tf: i for i, tf in enumerate(tfs)}

    # Check each entry in the mapping matrix against the annotated data
    for tf, regulation in tf_genes.items():
        if tf not in tf_idx:
            continue

        for gene, annotated in regulation.items():
            if gene not in gene_idx:
                continue

            predicted = A[gene_idx[gene], tf_idx[tf]]
            if annotated == 0:
                ambiguous += 1
            elif annotated > 0:
                total_pos += 1
                if predicted > 0:
                    correct_pos += 1
                elif predicted == 0:
                    dropped_pos += 1
            else:
                total_neg += 1
                if predicted < 0:
                    correct_neg += 1
                elif predicted == 0:
                    dropped_neg += 1

    # Print statistics
    annotated_neg = total_neg - dropped_neg
    annotated_pos = total_pos - dropped_pos
    correct = correct_neg + correct_pos
    dropped = dropped_neg + dropped_pos
    matched = annotated_neg + annotated_pos
    total = matched + dropped + ambiguous
    percent_dropped = 0 if total == 0 else 100*dropped/total
    percent_match = 0 if matched == 0 else 100*correct/matched
    percent_match_neg = 0 if annotated_neg == 0 else 100*correct_neg/annotated_neg
    percent_match_pos = 0 if annotated_pos == 0 else 100*correct_pos/annotated_pos
    print(f'Dropped regulation: {dropped}/{total} ({percent_dropped:.1f}%)')
    print(f'Overall matches: {correct}/{matched} ({percent_match:.1f}%)')
    print(f'Negative regulation matches: {correct_neg}/{annotated_neg} ({percent_match_neg:.1f}%)')
    print(f'Positive regulation matches: {correct_pos}/{annotated_pos} ({percent_match_pos:.1f}%)')
    print(f'Number of ambiguous regulatory interactions: {ambiguous}')

    # E prediction from results
    E_est = A.dot(P)
    no_prediction = E_est == 0
    n_zeros = np.sum(no_prediction)
    n_negatives = np.sum(E_est < 0)
    n_samples = E.shape[0] * E.shape[1]
    error = np.sqrt(np.mean(((E - E_est) / E)**2))
    error_predictions = np.sqrt(np.mean(((E[~no_prediction] - E_est[~no_prediction]) / E[~no_prediction])**2))
    print(f'\nError fitting original data: {error:.3f}')
    print(f'Error fitting original data (excluding unpredicted samples): {error_predictions:.3f}')
    print(f'No prediction made: {n_zeros}/{n_samples} ({100 * n_zeros / n_samples:.1f}%)')
    print(f'Negative prediction made: {n_negatives}/{n_samples} ({100 * n_negatives / n_samples:.1f}%)')

def plot_results(
        tf_genes: Dict[str, Dict[str, int]],
        A: np.ndarray,
        P: np.ndarray,
        genes: np.ndarray,
        tfs: np.ndarray,
        output_dir: str,
        ) -> None:
    """
    Plot NCA results for easier inspection.

    Args:
        tf_genes: relationship between TF and genes {TF: {gene: regulatory direction}}
        A: NCA solution for TF/gene matrix relationship (n genes, m TFs)
        P: NCA solution for TF/condition relationship (m TFs, o conditions)
        genes: IDs for each gene corresponding to rows in A (n genes)
        tfs: names of each TF corresponding to columns in A (m TFs)
        output_dir: path to directory to save the plot
    """

    def plot(ax, series, label, hist_range, n_bins, color):
        mean = sum(series) / len(series) if len(series) else 0
        ax.hist(series, color=color, bins=n_bins, range=hist_range, alpha=0.5, label=label)
        ax.axvline(mean, color=color, linestyle='--', label=f'{label} mean: {mean:.2f}')

    annotated_neg = []
    annotated_pos = []
    annotated_amb = []

    # Check each entry in the mapping matrix against the annotated data
    for i, j in zip(*np.where(A)):
        gene = genes[i]
        tf = tfs[j]
        annotated = tf_genes.get(tf, {}).get(gene)

        if annotated is not None:
            predicted = A[i, j]

            if annotated > 0:
                annotated_pos.append(predicted)
            elif annotated < 0:
                annotated_neg.append(predicted)
            else:
                annotated_amb.append(predicted)

    # Range for A histogram
    min_val = min(
        min(annotated_neg) if annotated_neg else 0,
        min(annotated_pos) if annotated_pos else 0,
        min(annotated_amb) if annotated_amb else 0,
        )
    max_val = max(
        max(annotated_neg) if annotated_neg else 0,
        max(annotated_pos) if annotated_pos else 0,
        max(annotated_amb) if annotated_amb else 0,
        )
    if min_val == max_val:
        max_val += 1
    hist_range = (
        np.floor(min_val),
        np.ceil(max_val),
    )
    n_bins = int(np.ceil(5*(hist_range[1] - hist_range[0])))
    cmap = plt.get_cmap('tab10')

    # TF regulation types
    positive_tfs = np.array([np.all([v > 0 for v in tf_genes.get(tf, {}).values()]) for tf in tfs])
    negative_tfs = np.array([np.all([v < 0 for v in tf_genes.get(tf, {}).values()]) for tf in tfs])
    combined_tfs = (~positive_tfs) & (~negative_tfs)

    # Plot figure
    plt.figure(figsize=(10, 20))
    gs = gridspec.GridSpec(4, 1)
    ax_A = plt.subplot(gs[0, 0])
    ax_P1 = plt.subplot(gs[1, 0])
    ax_P2 = plt.subplot(gs[2, 0])
    ax_P3 = plt.subplot(gs[3, 0])

    ## Plot A results
    plot(ax_A, annotated_neg, 'Negative', hist_range, n_bins, cmap(0))
    plot(ax_A, annotated_pos, 'Positive', hist_range, n_bins, cmap(1))
    plot(ax_A, annotated_amb, 'Ambiguous', hist_range, n_bins, cmap(2))
    ax_A.legend(fontsize=8, frameon=False)
    ax_A.set_xlabel('TF-Gene Interation\n(A entries)', fontsize=10)
    ax_A.set_ylabel('Count', fontsize=10)

    ## Plot P results
    P_ave = P.mean(1)
    hist_range = (np.floor(P_ave.min()), np.ceil(P_ave.max()))
    n_bins = 5 * int(hist_range[1] - hist_range[0])
    plot(ax_P1, P_ave[negative_tfs], 'All negative', hist_range, n_bins, cmap(0))
    plot(ax_P1, P_ave[positive_tfs], 'All positive', hist_range, n_bins, cmap(1))
    plot(ax_P1, P_ave[combined_tfs], 'Multiple', hist_range, n_bins, cmap(2))
    ax_P1.legend(fontsize=8, frameon=False)
    ax_P1.set_xlabel('TF Average Activity\n(mean of P rows)', fontsize=10)
    ax_P1.set_ylabel('Count', fontsize=10)

    ## Plot P range results
    P_range = P.max(1) - P.min(1)
    hist_range = (np.floor(P_range.min()), np.ceil(P_range.max()))
    n_bins = 2 * int(hist_range[1] - hist_range[0])
    plot(ax_P2, P_range[negative_tfs], 'All negative', hist_range, n_bins, cmap(0))
    plot(ax_P2, P_range[positive_tfs], 'All positive', hist_range, n_bins, cmap(1))
    plot(ax_P2, P_range[combined_tfs], 'Multiple', hist_range, n_bins, cmap(2))
    ax_P2.legend(fontsize=8, frameon=False)
    ax_P2.set_xlabel('TF Activity Range\n(range of P rows)', fontsize=10)
    ax_P2.set_ylabel('Count', fontsize=10)

    ## Plot P distributions
    for tf_activity in P:
        ax_P3.hist(tf_activity, linewidth=1, alpha=0.2, histtype='step')
    ax_P3.set_xlabel('TF Activity\n(P rows)', fontsize=10)
    ax_P3.set_ylabel('Count', fontsize=10)

    ## Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, HISTOGRAM_FILE))
    plt.close('all')

def parse_args() -> argparse.Namespace:
    """Parse command line args for options to run."""

    parser = argparse.ArgumentParser()

    default_nca = nca.METHODS[0]
    default_label = 'nca-results'
    default_iterative_iterations = 100
    default_iterative_split = 5
    default_robust_iterations = 100
    default_status_update = 0.1

    # General options
    parser.add_argument('-l', '--label',
        default=default_label,
        help=f'Label for output directory to save results to (default: {default_label}).')
    parser.add_argument('-v', '--verbose',
        action='store_true',
        help='If set, prints status updates for creating initial network map.')

    # Options for efficient analysis by reusing saved results
    parser.add_argument('-a', '--analysis',
        help='Path to directory containing saved regulation data to just run analysis on. Will skip NCA.')
    parser.add_argument('-c', '--cache',
        help='Path to cache directory to load network files. Defaults to output directory if not specified.')
    parser.add_argument('-f', '--force',
        action='store_true',
        help='If set, force a rerun of identifying the initial TF map, otherwise use cached values.')

    # Data options
    parser.add_argument('--linear',
        action='store_true',
        help='If set, use linear counts from sequencing data other keep log2 counts.')
    parser.add_argument('--global-expression',
        action='store_true',
        help='If set, creates pseudo transcription factors to capture global expression.')
    parser.add_argument('--split',
        action='store_true',
        help='If set, split transcription factors into positive and negative regulation.')

    # NCA options
    parser.add_argument('-m', '--method',
        choices=nca.METHODS,
        default=default_nca,
        help=f'NCA method to use, defined in nca.py (default: {default_nca}).')
    parser.add_argument('--status',
        type=float,
        default=default_status_update,
        help=f'Print status update every time this fraction of steps has been completed (default: {default_status_update}).')

    ## ISNCA specific options
    parser.add_argument('-i', '--iterative',
        action='store_true',
        help='If set, performs iterative sub-network component analysis.')
    parser.add_argument('--iterative-iterations',
        type=int,
        default=default_iterative_iterations,
        help=f'Maximum iterations for iterative sub-network component analysis (default: {default_iterative_iterations}).')
    parser.add_argument('--iterative-splits',
        type=int,
        default=default_iterative_split,
        help=f'Maximum splits for iterative sub-network component analysis (default: {default_iterative_split}).')

    ## ROBNCA specific options
    parser.add_argument('--robust-iterations',
        type=int,
        default=default_robust_iterations,
        help=f'Maximum iterations for robust_nca (default: {default_robust_iterations}).')

    return parser.parse_args()


if __name__ == '__main__':
    start = time.time()
    args = parse_args()

    # TODO: normalize seq data or only use EcoMAC data for now
    # TODO: handle options better when loading analysis or cache
    #   - seq_data can be different with args.linear
    #   - tf_genes can be different with args.split
    #   - tfs might not match saved regulation with args.global_expression
    print('Loading data from files...')
    seq_data = load_seq_data(args.linear)
    genes = load_gene_names()
    tf_genes = load_tf_gene_interactions(split=args.split, verbose=args.verbose)

    if args.analysis is None:
        # Check input cached files
        output_dir = os.path.join(OUTPUT_DIR, args.label)
        os.makedirs(output_dir, exist_ok=True)
        cache_dir = args.cache if args.cache else os.path.join(output_dir, 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        network_cache_file = os.path.join(cache_dir, NETWORK_CACHE_FILE)
        tf_cache_file = os.path.join(cache_dir, TF_CACHE_FILE)
        no_cache = not (os.path.exists(network_cache_file) and os.path.exists(tf_cache_file))

        # Create or load network mapping and TF IDs
        if args.force or no_cache:
            print('Creating initial network mapping...')
            initial_tf_map, tfs = create_tf_map(genes, tf_genes, verbose=args.verbose)

            if not args.iterative:
                initial_tf_map, tfs = nca.nca_criteria_check(initial_tf_map, tfs, verbose=args.verbose)
        else:
            print('Loading cached initial network mapping...')
            initial_tf_map = np.load(network_cache_file)
            tfs = np.load(tf_cache_file)

        # Output cached files
        cache_dir = os.path.join(output_dir, 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        network_cache_file = os.path.join(cache_dir, NETWORK_CACHE_FILE)
        tf_cache_file = os.path.join(cache_dir, TF_CACHE_FILE)
        np.save(network_cache_file, initial_tf_map)
        np.save(tf_cache_file, tfs)

        if args.global_expression:
            tfs, initial_tf_map = add_global_expression(tfs, mapping=initial_tf_map)

        # Solve NCA problem
        nca_method = getattr(nca, args.method)
        if args.iterative:
            A, P, tfs = nca.iterative_sub_nca(nca_method, seq_data, initial_tf_map, tfs,
                n_iters=args.iterative_iterations, splits=args.iterative_splits,
                robust_iters=args.robust_iterations, status_step=args.status, verbose=args.verbose)
        else:
            A, P = nca_method(seq_data, initial_tf_map, n_iters=args.robust_iterations, status_step=args.status)

        # Save results
        save_regulation(A, P, genes, tfs, output_dir)
    else:
        print('Loading regulation results without running NCA...')
        output_dir = args.analysis
        if not os.path.exists(output_dir):
            raise IOError(f'Directory does not exist: {output_dir}')
        cache_dir = args.cache if args.cache else os.path.join(args.analysis, 'cache')

        tfs = np.load(os.path.join(cache_dir, TF_CACHE_FILE))
        if args.global_expression:
            tfs, _ = add_global_expression(tfs)

        A, P = load_regulation(args.analysis, genes, tfs)

    # Assess results of analysis
    match_statistics(seq_data, A, P, tf_genes, genes, tfs)
    plot_results(tf_genes, A, P, genes, tfs, output_dir)

    print(f'Completed in {(time.time() - start) / 60:.1f} min')
