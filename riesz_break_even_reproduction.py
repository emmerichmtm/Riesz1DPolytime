
"""
Reproducibility benchmark for fixed-cardinality 1-D Riesz s-energy subset selection.

Compares:
1. Complete enumeration over all k-subsets.
2. Explicit threshold min-cut construction from the paper.

Instances:
    n = 2k
    x_i = i + 0.05 sin(1.7 i) + 0.01 (i/n)^2, i=0,...,n-1
    s = 1

Outputs:
    riesz_break_even_benchmark.csv
    riesz_break_even_times.png

Requirements:
    numpy pandas matplotlib networkx numba
"""

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from numba import njit


@njit(cache=False)
def exhaustive_numba(W, n, k):
    """Complete enumeration using a lexicographic combination generator."""
    comb = np.empty(k, dtype=np.int64)
    best = np.empty(k, dtype=np.int64)
    for i in range(k):
        comb[i] = i
        best[i] = i

    best_val = 1.0e300
    count = 0

    while True:
        e = 0.0
        for a in range(k):
            ia = comb[a]
            for b in range(a + 1, k):
                e += W[ia, comb[b]]

        if e < best_val:
            best_val = e
            for i in range(k):
                best[i] = comb[i]

        count += 1

        idx = k - 1
        while idx >= 0 and comb[idx] == idx + n - k:
            idx -= 1

        if idx < 0:
            break

        comb[idx] += 1
        for j in range(idx + 1, k):
            comb[j] = comb[j - 1] + 1

    return best_val, best, count


def make_points(n):
    """Deterministic mildly irregular ordered point set."""
    i = np.arange(n, dtype=float)
    return i + 0.05 * np.sin(1.7 * i) + 0.01 * (i / n) ** 2


def pairwise_weights(x, s=1.0):
    n = len(x)
    W = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            W[i, j] = 1.0 / ((x[j] - x[i]) ** s)
    return W


def energy_subset(x, idx, s=1.0):
    e = 0.0
    idx = list(idx)
    for a, i in enumerate(idx):
        for j in idx[a + 1:]:
            e += (x[j] - x[i]) ** (-s)
    return e


def mincut_riesz(x, k, s=1.0):
    """
    Explicit threshold min-cut construction.

    Variables:
        z_{r,t} = 1[y_r >= t], r=0,...,k-1, t=0,...,m-1
    Feasible subset:
        i_r = r + y_r
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    m = n - k

    if k <= 1:
        return 0.0, tuple(range(k)), 2, 0
    if k == n:
        return energy_subset(x, range(n), s), tuple(range(n)), 2, 0

    unary = np.zeros((k, m), dtype=float)
    pair_edges = {}

    def V(p, q, a, b):
        return (x[q + b] - x[p + a]) ** (-s)

    for p in range(k):
        for q in range(p + 1, k):
            beta = np.zeros(m + 1, dtype=float)

            for u in range(1, m + 1):
                beta[u] = V(p, q, 0, u) - V(p, q, 0, u - 1)
                unary[q, u - 1] += beta[u]

            Delta = np.zeros((m + 1, m + 1), dtype=float)
            for t in range(1, m + 1):
                for u in range(t + 1, m + 1):
                    Delta[t, u] = (
                        V(p, q, t, u)
                        - V(p, q, t - 1, u)
                        - V(p, q, t, u - 1)
                        + V(p, q, t - 1, u - 1)
                    )

            for a in range(1, m + 1):
                alpha = V(p, q, a, a) - V(p, q, a - 1, a - 1) - beta[a]
                for t in range(1, a):
                    alpha -= Delta[t, a]
                unary[p, a - 1] += alpha

            for t in range(1, m + 1):
                for u in range(t + 1, m + 1):
                    d = Delta[t, u]
                    if d < -1e-15:
                        w = -d
                        unary[p, t - 1] -= w
                        key = ((p, t - 1), (q, u - 1))
                        pair_edges[key] = pair_edges.get(key, 0.0) + w

    G = nx.DiGraph()
    S = "S"
    T = "T"
    threshold_nodes = [(r, t) for r in range(k) for t in range(m)]
    G.add_nodes_from(threshold_nodes + [S, T])

    finite_sum = 0.0

    for (a, b), w in pair_edges.items():
        w = float(w)
        G.add_edge(a, b, capacity=w)
        finite_sum += w

    for r in range(k):
        for t in range(m):
            c = float(unary[r, t])
            if c >= 0:
                if c > 1e-15:
                    G.add_edge((r, t), T, capacity=c)
                    finite_sum += c
            else:
                G.add_edge(S, (r, t), capacity=-c)
                finite_sum += -c

    U = finite_sum + 1.0

    for r in range(k):
        for t in range(m - 1):
            G.add_edge((r, t + 1), (r, t), capacity=U)

    for r in range(k - 1):
        for t in range(m):
            G.add_edge((r, t), (r + 1, t), capacity=U)

    _, (reachable, _) = nx.minimum_cut(
        G, S, T, capacity="capacity", flow_func=nx.algorithms.flow.preflow_push
    )

    y = []
    for r in range(k):
        y.append(sum((r, t) in reachable for t in range(m)))

    idx = tuple(r + y[r] for r in range(k))
    return energy_subset(x, idx, s), idx, G.number_of_nodes(), G.number_of_edges()


def run_benchmark():
    s = 1.0

    # Exhaustive enumeration is run until after practical break-even.
    exhaustive_ns = list(range(4, 32, 2))

    # Min-cut is run well beyond the break-even point.
    mincut_ns = list(range(4, 32, 2)) + [36, 40, 44, 48, 52, 56, 60]

    # Warm up Numba compilation outside timings.
    x0 = make_points(8)
    W0 = pairwise_weights(x0, s)
    exhaustive_numba(W0, 8, 4)

    rows = []

    for n in mincut_ns:
        k = n // 2
        x = make_points(n)
        W = pairwise_weights(x, s)
        comb_count = math.comb(n, k)

        brute_time = None
        brute_energy = None
        brute_subset = None
        if n in exhaustive_ns:
            t0 = time.perf_counter()
            brute_energy, brute_subset_arr, count = exhaustive_numba(W, n, k)
            brute_time = time.perf_counter() - t0
            brute_subset = tuple(int(i) for i in brute_subset_arr)
            assert count == comb_count

        t0 = time.perf_counter()
        cut_energy, cut_subset, nodes, arcs = mincut_riesz(x, k, s)
        cut_time = time.perf_counter() - t0

        if brute_energy is not None:
            assert abs(brute_energy - cut_energy) < 1e-8

        rows.append({
            "n": n,
            "k": k,
            "m": n - k,
            "comb_count": comb_count,
            "enum_time_s": brute_time,
            "mincut_time_s": cut_time,
            "speedup_enum_over_mincut": None if brute_time is None else brute_time / cut_time,
            "mincut_nodes": nodes,
            "mincut_arcs": arcs,
            "energy": cut_energy,
            "subset_zero_based": str(cut_subset),
            "subset_one_based": str(tuple(i + 1 for i in cut_subset)),
            "enum_energy": brute_energy,
            "enum_subset_zero_based": None if brute_subset is None else str(brute_subset),
            "theory_enum_work_k2C": (k * k) * comb_count,
            "theory_cut_balanced_k8": k ** 8,
            "theory_cut_user_n8_over_2": (n ** 8) / 2.0,
        })

    df = pd.DataFrame(rows)

    df.to_csv("riesz_break_even_benchmark.csv", index=False)

    measured = df[df["enum_time_s"].notna()].copy()

    plt.figure(figsize=(8, 5))
    plt.plot(measured["n"], measured["enum_time_s"], marker="o", label="complete enumeration")
    plt.plot(df["n"], df["mincut_time_s"], marker="s", label="threshold min-cut")
    plt.yscale("log")
    plt.xlabel("n with k = n/2")
    plt.ylabel("CPU time in seconds, log scale")
    plt.title("Riesz subset selection: enumeration vs. min-cut")
    plt.legend()
    plt.tight_layout()
    plt.savefig("riesz_break_even_times.png", dpi=200)

    return df


if __name__ == "__main__":
    df = run_benchmark()
    print(df.to_string(index=False))
