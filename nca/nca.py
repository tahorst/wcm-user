#! /usr/bin/env python

"""
NCA methods to use to solve E = AP given E and specified network connections in A.
"""

import numpy as np
import scipy.linalg


# Function names of the different NCA methods that have been implemented below
METHODS = ['robust_nca', 'fast_nca']


def nca_criteria_check(A: np.ndarray, tfs: np.ndarray, verbose: bool = True) -> (np.ndarray, np.ndarray):
    """
    Check criteria for the A matrix of NCA (E = AP).
    - full column rank in A
    - removing column and associated rows leads to rank of m-1 for m columns in A

    Args:
        A: initial A matrix for NCA algorithm (n genes, m TFs)
        tfs: IDs of transcription factors for each column in A
        verbose: If set, prints updates about removed columns

    Returns:
        A: updated A matrix with columns removed to satisfy NCA criteria
        tfs: updated TF IDs corresponding to the new A
    """

    if verbose:
        print(f'A shape: {A.shape}')

    # Filter out TFs that did not match any known genes for regulation (empty columns)
    col_sum = np.sum(A != 0, axis=0)
    if np.any(col_sum == 0):
        tf_idx_with_regulation = np.where(col_sum != 0)[0]
        if verbose:
            print('Removing empty columns for TFs:')
            for tf in tfs[col_sum == 0]:
                print(f'\t{tf}')
        A, tfs = nca_criteria_check(A[:, tf_idx_with_regulation], tfs[tf_idx_with_regulation], verbose=verbose)

    n_cols = A.shape[1]

    # Check if not full rank A matrix
    rank = np.linalg.matrix_rank(A)
    if rank < n_cols:
        if verbose:
            print('A matrix is not full rank because of TFs:')
        dependent = []
        for i in range(n_cols):
            mask = np.ones(n_cols, bool)
            mask[i] = False
            sub_rank = np.linalg.matrix_rank(A[:, mask])
            if sub_rank == rank:
                dependent.append(i)
                if verbose:
                    print(f'\t{tfs[i]}')

        # Check for dependent columns with LU decomposition and remove
        # _, _, U = scipy.linalg.lu(A)
        # dependent = np.array([i for i in range(A.shape[1]) if not np.any(U[i, :])])
        # TODO: better way of checking this or combine dependent columns into one to not lose all info
        dependent = np.array(dependent)
        mask = np.ones(n_cols, bool)
        mask[dependent] = False
        A, tfs = nca_criteria_check(A[:, mask], tfs[mask], verbose=verbose)
        rank = np.linalg.matrix_rank(A)
        n_cols = A.shape[1]

    # Check that reduced matrices have expected rank
    for i in range(n_cols):
        row_mask = A[:, i] == 0
        col_mask = np.ones(n_cols, bool)
        col_mask[i] = False

        reduced_A = A[row_mask, :][:, col_mask]
        sub_rank = np.linalg.matrix_rank(reduced_A)
        if sub_rank != rank - 1:
            # Filter out empty columns with rows removed
            col_sum = np.sum(A[row_mask, :] != 0, axis=0)
            tf_idx_with_regulation = col_sum != 0
            tf_idx_with_regulation[i] = True
            if not np.all(tf_idx_with_regulation):
                # TODO: turn into function with same code above
                if verbose:
                    print('Removing empty columns from reduced matrix for TFs:')
                    for tf in tfs[~tf_idx_with_regulation]:
                        print(f'\t{tf}')
                A, tfs = nca_criteria_check(A[:, tf_idx_with_regulation], tfs[tf_idx_with_regulation], verbose=verbose)
            else:
                # TODO; turn into function with same code above
                if verbose:
                    print('Reduced matrix is not full rank because of TFs:')
                rows_removed = A[row_mask, :]
                dependent = []
                for j in range(n_cols):
                    col_mask = np.ones(n_cols, bool)
                    col_mask[j] = False
                    new_rank = np.linalg.matrix_rank(rows_removed[:, col_mask])
                    if new_rank == sub_rank:
                        dependent.append(j)
                        if verbose:
                            print(f'\t{tfs[j]}')

                dependent = np.array(dependent)
                mask = np.ones(n_cols, bool)
                mask[dependent] = False
                A, tfs = nca_criteria_check(A[:, mask], tfs[mask], verbose=verbose)
            break

    return A, tfs

def fast_nca(E: np.ndarray, A: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Perform FastNCA on dataset E with network connectivity specified by A for the
    problem: E = AP. Based on matlab implementation from Chang et al. 2008.
    http://www.eee.hku.hk/~cqchang/FastNCA.m

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    print('Solving with FastNCA...')
    n_cols = A.shape[1]
    U, S, V = scipy.linalg.svd(E)
    U = U[:, :n_cols]
    S = np.diag(S[:n_cols])
    V = V[:, :n_cols].T
    E_approx = U.dot(S).dot(V)

    A_hat = np.zeros_like(A)
    status_step = 0.1
    next_update = status_step
    for i in range(n_cols):
        if (i + 1) / n_cols >= next_update:
            complete = np.floor((i + 1) / n_cols * 100)
            print(f'{complete:.0f}% complete...')
            next_update += status_step

        U0 = U[A[:, i] == 0, :]
        if U0.shape[0] < n_cols:
            _, _, M = scipy.linalg.svd(U0)
        else:
            M, _, _ = scipy.linalg.svd(U0.T)
        v = M[:, -1]
        a = U.dot(v)
        a[A[:, i] == 0] = 0
        A_hat[:, i] = a / np.sum(np.abs(a)) * np.sum(A[:, i] != 0)

    P_hat = np.linalg.lstsq(A_hat, E_approx, rcond=None)[0]

    return A_hat, P_hat

def robust_nca(E: np.ndarray, A: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Perform ROBNCA on dataset E with network connectivity specified by A for the
    problem: E = AP. Based on method in Noor et al. Bioinformatics. 2013.

    Args:
        E: data to solve NCA for (n genes, m conditions)
        A: network connectivity (n genes, o TFs)

    Returns:
        A_est: estimated A based fit to data (n genes, o TFs)
        P_est: estimated P based fit to data (o TFs, m conditions)
    """

    print('Solving with ROBNCA...')
    n_tfs = A.shape[1]
    A_est = A.copy()
    outliers = np.zeros_like(E)
    zero_mask = A == 0

    lambda_ = 0.7  # parameter that can be tuned for sparsity to determine outliers
    n_iters = 100
    status_step = 0.1
    next_update = status_step
    for it in range(n_iters):
        if (it + 1) / n_iters >= next_update:
            complete = np.floor((it + 1) / n_iters * 100)
            print(f'{complete:.0f}% complete...')
            next_update += status_step

        # Update S
        X = E - outliers
        S = np.linalg.inv(A_est.T.dot(A_est)).dot(A_est.T).dot(X)

        # Update A
        Qinv = np.linalg.inv(S.dot(S.T))
        for col, a in enumerate(A_est):
            zero_idx = np.where(a == 0)[0]
            Cn = np.zeros((len(zero_idx), n_tfs))
            Cn[range(len(zero_idx)), zero_idx] = 1
            psi = Cn.dot(Qinv).dot(Cn.T)
            w = S.dot(X.T[:, col])
            A_est[col, :] = Qinv.dot(w - Cn.T.dot(np.linalg.inv(psi)).dot(Cn).dot(Qinv).dot(w))
        A_est[zero_mask] = 0

        # Update outliers
        B = E - A_est.dot(S)
        norm = np.linalg.norm(B, axis=0)
        outliers = B * np.fmax(0, (norm - lambda_ / 2)) / norm

    P_est = S
    return A_est, P_est
