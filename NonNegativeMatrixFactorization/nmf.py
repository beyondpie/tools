"""
Non-negative matrix factorization (NMF)

Original Author: Yang Eric Li
Updated by: Songpeng Zu
Updated time: 2024-08-23

- remove the unused inputs for xgi and ygi
- remove the unused method lda
- Add another input support:
  h5 file (mat, feature by cluster)
"""

import argparse
import h5py
import numpy as np
from scipy import sparse
from scipy.sparse import load_npz
from sklearn.decomposition import NMF
from time import perf_counter as pc
import matplotlib.pyplot as plt

plt.switch_backend("agg")

parser = argparse.ArgumentParser(description="Run NMF using sklearn.")
parser.add_argument(
    "-i",
    "--inputF",
    type=str,
    dest="inputF",
    help="input matrix in npz format",
    default=None,
)

parser.add_argument("-r", "--rank", type=int, dest="rank", help="an integer for rank ")
parser.add_argument(
    "-n", "--seed", type=int, dest="seed", default=1)
parser.add_argument(
    "-o", "--outPrefix", type=str, dest="outPrefix", help="output prefix"
)
parser.add_argument("-p", "--pbych5", type=str, default=None)

args = parser.parse_args()


def read_npz(inputf):
    """
    Read snATAC data in npz format. The matrix's shape is ### (bins) x ### (samples).
    It contains only positive data (boolean).

    Return the dense data matrix.
    """
    V = load_npz(inputf)
    V = V.tocsr()
    return V


def saveH(prefix, X):
    print("=== write matrix H ===")
    fileN = [prefix, "H", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def saveW(prefix, X):
    print("=== write matrix W ===")
    fileN = [prefix, "W", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def run_nmf_single(V, rank, n, prefix):
    """
    Run standard NMF on data set.

    :param V: Target matrix with gene expression data.
    :type V: `sparse.matrix`
    :param rank: Factorization rank.
    :type rank: `int`
    :param n: randome state
    """
    print("perform NMF 2 in nndsvd model")
    model = NMF(
        n_components=rank,
        init="nndsvd",
        random_state=n,
        verbose=True,
        max_iter=500,
    )
    W = model.fit_transform(V)
    H = model.components_
    saveH(prefix, H)
    saveW(prefix, W)

if __name__ == "__main__":
    """Run standard NMF on rank"""
    outPrefix = args.outPrefix
    rank = args.rank
    n = args.seed
    start_time = pc()
    if args.inputF is not None:
        print(f"Loading mat from npz file {args.inputF}.")
        V = read_npz(args.inputF)
    elif args.pbych5 is not None:
        print(f"Loading mat from h5 file {args.pbych5}.")
        with h5py.File(args.pbych5, mode="r") as f:
            t = f["X"]["mat"][()]
        V = sparse.csr_matrix(t)
        del t
    else:
        raise FileNotFoundError("Lack of input mat.")
    run_nmf_single(V, rank, n, outPrefix)
    end_time = pc()
    print("Used (secs): ", end_time - start_time)
