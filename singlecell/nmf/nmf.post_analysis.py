"""
After running nmf independently, we read the results and summize them.

Original Author: Yang Eric Li.
Updated by: Songpeng Zu
Updated time: 2024-08-23
"""

from warnings import warn
from os.path import exists, basename
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from time import perf_counter as pc
import matplotlib.pyplot as plt
import argparse
import itertools
import math

parser = argparse.ArgumentParser(description="Run NMF using sklearn.")
parser.add_argument("-i", "--inputPrefix", type=str, dest="input_prefix")
parser.add_argument("-r", "--rank", type=int, dest="rank")
parser.add_argument("-d", "--defaultRandomState", type=int, dest="default_random_state")
parser.add_argument("-s", "--startRandomState", type=int, dest="start_random_state")
parser.add_argument("-e", "--endRandomState", type=int, dest="end_random_state")
parser.add_argument("-o", "--outPrefix", type=str, dest="out_prefix")

args = parser.parse_args()

plt.switch_backend("agg")
try:
    from matplotlib.pyplot import imshow, set_cmap
except ImportError:
    warn("Matplotlib must be installed.")


def saveC(prefix, X):
    print("=== write matrix C ===")
    fileN = [prefix, "C", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def cal_featureScore_kim(W):
    """extract feature from W"""
    print("=== extract feature from W ===")
    k = W.shape[1]
    m = W.shape[0]
    s_list = []
    for i in range(m):
        rowsum = np.sum(W[i,])
        p_iq_x_list = []
        for q in range(k):
            p_iq = W[i, q] / rowsum
            if p_iq != 0:
                tmp = p_iq * math.log(p_iq, 2)
            else:
                tmp = 0
            p_iq_x_list.append(tmp)
        s = 1 + 1 / math.log(k, 2) * np.sum(p_iq_x_list)
        s_list.append(s)
    return s_list


def predict_H(H):
    """extract feature from H"""
    print("=== extract feature from H ===")
    colmax = np.amax(H, axis=0)
    colsum = np.sum(H, axis=0)
    p = colmax / colsum
    idx = H.argmax(axis=0)
    out = [idx, p]
    return out


def cal_connectivity(H, idx):
    """calculate connectivity matrix"""
    print("=== calculate connectivity matrix ===")
    connectivity_mat = np.zeros((H.shape[1], H.shape[1]))
    classN = H.shape[0]
    for i in range(classN):
        xidx = list(np.concatenate(np.where(idx == i)))
        iterables = [xidx, xidx]
        for t in itertools.product(*iterables):
            connectivity_mat[t[0], t[1]] = 1
    return connectivity_mat


def post_nmf(
    input_prefix,
    n_mods,
    out_prefix,
    start_random_state=0,
    end_random_state=9,
    default_random_state=0,
):
    # load all the generated H
    H_candidates = [
        f"{input_prefix}.r{n_mods}.n{i}.H.mx"
        for i in range(start_random_state, end_random_state)
    ]
    H_exists = [j for j in H_candidates if exists(j)]
    print("{len(H_exists)} H.mx files are found for r{n_mods}.")
    used_random_states = [basename(i).split(".")[3] for i in H_exists]
    print(f"Found random states: {','.join(used_random_states)}")
    # load default H and W
    H = np.loadtxt(f"{input_prefix}.r{n_mods}.n{default_random_state}.H.mx")
    W = np.loadtxt(f"{input_prefix}.r{n_mods}.n{default_random_state}.W.mx")
    # get consensus matrix (depend on Hs)
    n_cluster = H.shape[1]
    consensus = np.zeros(shape=(n_cluster, n_cluster))
    for ih in H_exists:
        iH = np.loadtxt(ih)
        consensus += cal_connectivity(iH, predict_H(iH)[0])
    consensus /= len(H_exists)
    # plot and save results
    p_consensus = reorder(consensus)
    outer = f"{out_prefix}.r{n_mods}"
    plotC(outer, p_consensus, n_mods)
    saveC(outer, p_consensus)
    o_fsW = cal_featureScore_kim(W)
    o_predH = predict_H(H)
    np.savetxt(".".join([outer, "featureScore_W.txt"]), o_fsW)
    np.savetxt(
        ".".join([outer, "predict_H.txt"]), np.squeeze(o_predH).T, delimiter="\t"
    )


def plotC(prefix, C, rank):
    """
    Plot reordered consensus matrix.

    :param C: Reordered consensus matrix.
    :type C: numpy.ndarray`
    :param rank: Factorization rank.
    :type rank: `int`
    """
    fig = plt.figure(figsize=(5, 5), dpi=100)
    imshow(C)
    set_cmap("RdBu_r")
    fileN = [prefix, "C", "png"]
    fileN = ".".join(fileN)
    fig.savefig(fileN)


def reorder(C):
    """
    Reorder consensus matrix.

    :param C: Consensus matrix.
    :type C: `numpy.ndarray`
    """
    Y = 1 - C
    Z = linkage(squareform(Y), method="average")
    ivl = leaves_list(Z)
    ivl = ivl[::-1]
    return C[:, ivl][ivl, :]


if __name__ == "__main__":
    """Run post analysis to summaize nmf."""
    input_prefix = args.input_prefix
    out_prefix = args.out_prefix
    rank = args.rank
    default_random_state = args.default_random_state
    start_random_state = args.start_random_state
    end_random_state = args.end_random_state
    start_time = pc()
    post_nmf(
        input_prefix=input_prefix,
        n_mods=rank,
        out_prefix=out_prefix,
        start_random_state=start_random_state,
        end_random_state=end_random_state,
        default_random_state=default_random_state,
    )
    end_time = pc()
    print("Used (secs): ", end_time - start_time)
