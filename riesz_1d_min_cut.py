#!/usr/bin/env python3
"""
Polynomial-time graph-cut algorithm for the one-dimensional fixed-cardinality
minimum Riesz s-energy subset problem.

Problem
-------
Given sorted real points x_0 < ... < x_{n-1}, a cardinality k, and s > 0,
choose indices i_0 < ... < i_{k-1} minimizing

    sum_{p<q} (x[i_q] - x[i_p])^{-s}.

The implementation uses the threshold-variable min-cut formulation described
in the accompanying README.  It has no third-party dependencies.

The numerical implementation uses floating point capacities.  For  exact arithmetic with rational points and fixed integer s, the same
combinatorial graph can be built with rational capacities and solved by an
exact max-flow implementation.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import combinations
import argparse
import math
import random
from typing import Iterable, List, Sequence, Tuple


@dataclass
class Result:
    indices: List[int]
    points: List[float]
    energy: float


class Dinic:
    """A small directed max-flow implementation for nonnegative capacities."""

    def __init__(self, n: int):
        self.n = n
        self.adj: List[List[List[float | int | list]]] = [[] for _ in range(n)]

    def add_edge(self, u: int, v: int, c: float) -> None:
        if c <= 1e-15:
            return
        fwd = [v, float(c), None]
        rev = [u, 0.0, fwd]
        fwd[2] = rev
        self.adj[u].append(fwd)
        self.adj[v].append(rev)

    def maxflow(self, source: int, sink: int) -> float:
        flow = 0.0
        n = self.n

        while True:
            level = [-1] * n
            level[source] = 0
            q = deque([source])
            while q:
                u = q.popleft()
                for e in self.adj[u]:
                    v, cap, _ = e
                    if cap > 1e-12 and level[v] < 0:
                        level[v] = level[u] + 1
                        q.append(v)
            if level[sink] < 0:
                break

            it = [0] * n

            def dfs(u: int, f: float) -> float:
                if u == sink:
                    return f
                while it[u] < len(self.adj[u]):
                    e = self.adj[u][it[u]]
                    v, cap, rev = e
                    if cap > 1e-12 and level[v] == level[u] + 1:
                        pushed = dfs(v, min(f, cap))
                        if pushed > 1e-12:
                            e[1] -= pushed
                            rev[1] += pushed
                            return pushed
                    it[u] += 1
                return 0.0

            while True:
                pushed = dfs(source, 1e100)
                if pushed <= 1e-12:
                    break
                flow += pushed

        return flow

    def source_side(self, source: int) -> List[bool]:
        """Return vertices reachable from source in the residual graph."""
        seen = [False] * self.n
        seen[source] = True
        stack = [source]
        while stack:
            u = stack.pop()
            for v, cap, _ in self.adj[u]:
                if cap > 1e-12 and not seen[v]:
                    seen[v] = True
                    stack.append(v)
        return seen


def riesz_energy(points: Sequence[float], indices: Sequence[int], s: float) -> float:
    """Evaluate the Riesz s-energy of the selected indices."""
    return sum((points[j] - points[i]) ** (-s) for i, j in combinations(indices, 2))


def _pair_decomposition(
    x: Sequence[float], s: float, p: int, q: int, m: int
) -> Tuple[float, List[float], List[float], List[List[float]]]:
    """
    Decompose V(a,b)=(x[q+b]-x[p+a])^{-s}, 0<=a<=b<=m, into threshold terms.

    V(a,b) = C + sum_{t<=a} alpha[t] + sum_{u<=b} beta[u]
             + sum_{t<=a,u<=b,t<u} Delta[t][u].

    The Monge inequality implies Delta[t][u] <= 0 for t<u.
    """

    def V(a: int, b: int) -> float:
        return (x[q + b] - x[p + a]) ** (-s)

    C = V(0, 0)
    alpha = [0.0] * (m + 1)
    beta = [0.0] * (m + 1)
    delta = [[0.0] * (m + 1) for _ in range(m + 1)]

    for u in range(1, m + 1):
        beta[u] = V(0, u) - V(0, u - 1)

    for t in range(1, m + 1):
        for u in range(t + 1, m + 1):
            delta[t][u] = V(t, u) - V(t - 1, u) - V(t, u - 1) + V(t - 1, u - 1)
            if delta[t][u] > 1e-9:
                raise ArithmeticError(
                    f"Monge sign failed numerically for p={p}, q={q}, t={t}, u={u}: "
                    f"Delta={delta[t][u]}"
                )

    # Diagonal recurrence.  Since V(a,a)-V(a-1,a-1) = alpha[a]+beta[a]
    # + sum_{t=1}^{a-1} Delta[t][a], this determines alpha[a].
    for a in range(1, m + 1):
        diagonal_delta = sum(delta[t][a] for t in range(1, a))
        alpha[a] = V(a, a) - V(a - 1, a - 1) - beta[a] - diagonal_delta

    return C, alpha, beta, delta


def minimize_riesz_1d(points: Iterable[float], k: int, s: float = 1.0) -> Result:
    """
    Minimize one-dimensional fixed-cardinality Riesz s-energy by one min-cut.

    Parameters
    ----------
    points:
        Iterable of distinct real coordinates.  They are sorted internally.
    k:
        Required cardinality.
    s:
        Positive Riesz exponent.

    Returns
    -------
    Result(indices, points, energy), where indices are zero-based indices in the
    internally sorted point list.
    """
    if s <= 0:
        raise ValueError("s must be positive")

    x = sorted(float(v) for v in points)
    n = len(x)
    if len(set(x)) != n:
        raise ValueError("points must be distinct")
    if not (0 <= k <= n):
        raise ValueError("need 0 <= k <= n")

    if k == 0:
        return Result([], [], 0.0)
    if k == 1:
        return Result([0], [x[0]], 0.0)
    if k == n:
        idx = list(range(n))
        return Result(idx, [x[i] for i in idx], riesz_energy(x, idx, s))

    m = n - k
    num_nodes = k * m
    source = num_nodes
    sink = num_nodes + 1
    graph = Dinic(num_nodes + 2)
    unary = [0.0] * num_nodes

    def node(r: int, t: int) -> int:
        # r is zero-based rank; t is a threshold in {1,...,m}.
        return r * m + (t - 1)

    min_gap = min(x[i + 1] - x[i] for i in range(n - 1))
    upper_energy = math.comb(k, 2) * min_gap ** (-s)
    inf_capacity = max(1e12, 1000.0 * upper_energy + 1.0)

    # Closure constraints z_{r,t+1} <= z_{r,t}.
    for r in range(k):
        for t in range(1, m):
            graph.add_edge(node(r, t + 1), node(r, t), inf_capacity)

    # Monotonicity constraints y_r <= y_{r+1}, i.e. z_{r,t} <= z_{r+1,t}.
    for r in range(k - 1):
        for t in range(1, m + 1):
            graph.add_edge(node(r, t), node(r + 1, t), inf_capacity)

    # Pairwise Riesz terms.
    for p in range(k):
        for q in range(p + 1, k):
            _, alpha, beta, delta = _pair_decomposition(x, s, p, q, m)
            for t in range(1, m + 1):
                unary[node(p, t)] += alpha[t]
                unary[node(q, t)] += beta[t]
            for t in range(1, m + 1):
                for u in range(t + 1, m + 1):
                    d = delta[t][u]
                    if d < -1e-14:
                        w = -d
                        # -w xy = -w x + w x(1-y).
                        unary[node(p, t)] -= w
                        graph.add_edge(node(p, t), node(q, u), w)

    # Unary source/sink terms.
    constant = 0.0
    for v, coeff in enumerate(unary):
        if coeff >= 0:
            graph.add_edge(v, sink, coeff)
        else:
            constant += coeff
            graph.add_edge(source, v, -coeff)

    graph.maxflow(source, sink)
    side = graph.source_side(source)

    y = []
    for r in range(k):
        yr = sum(1 for t in range(1, m + 1) if side[node(r, t)])
        y.append(yr)

    indices = [r + y[r] for r in range(k)]
    selected = [x[i] for i in indices]
    return Result(indices, selected, riesz_energy(x, indices, s))


def brute_force_minimize_riesz_1d(points: Iterable[float], k: int, s: float = 1.0) -> Result:
    """Exhaustive checker for small instances."""
    x = sorted(float(v) for v in points)
    if not (0 <= k <= len(x)):
        raise ValueError("need 0 <= k <= n")
    best_indices: List[int] | None = None
    best_energy = float("inf")
    for indices_tuple in combinations(range(len(x)), k):
        energy = riesz_energy(x, indices_tuple, s)
        if energy < best_energy:
            best_energy = energy
            best_indices = list(indices_tuple)
    assert best_indices is not None or k == 0
    if k == 0:
        return Result([], [], 0.0)
    return Result(best_indices, [x[i] for i in best_indices], best_energy)


EXAMPLES = [
    ([0, 1, 2, 3, 4], 3, 1.0),
    ([0, 0.4, 1.1, 2.8, 3.0, 5.0], 3, 2.0),
    ([0, 1, 10, 11, 12], 3, 1.0),
    ([0, 1, 2, 4, 7, 11], 4, 1.5),
    ([0, 1, 2, 3, 10, 11, 12, 20], 4, 1.0),
    ([0, 0.2, 0.9, 2.7, 4.1, 4.2, 8.0], 3, 3.0),
    ([0, 5, 6, 7, 8, 20, 21, 40], 4, 0.5),
    ([0, 1, 1.5, 2.2, 6.0, 9.0, 9.1, 14.0, 18.0], 5, 2.0),
]


def run_examples() -> None:
    for x, k, s in EXAMPLES:
        alg = minimize_riesz_1d(x, k, s)
        brute = brute_force_minimize_riesz_1d(x, k, s)
        ok = abs(alg.energy - brute.energy) <= 1e-8 * max(1.0, abs(brute.energy))
        one_based = [i + 1 for i in alg.indices]
        print(f"x={x}, k={k}, s={s}")
        print(f"  min-cut indices 0-based: {alg.indices}")
        print(f"  min-cut indices 1-based: {one_based}")
        print(f"  selected points: {alg.points}")
        print(f"  energy: {alg.energy:.12g}")
        print(f"  agrees with brute force: {ok}")
        print()


def self_test(max_n: int = 10, trials: int = 100, seed: int = 1) -> None:
    random.seed(seed)
    for n in range(2, max_n + 1):
        for k in range(1, n + 1):
            for _ in range(trials):
                x = sorted(random.sample(range(0, 200), n))
                s = random.choice([0.5, 1.0, 1.5, 2.0, 3.0])
                alg = minimize_riesz_1d(x, k, s)
                brute = brute_force_minimize_riesz_1d(x, k, s)
                if abs(alg.energy - brute.energy) > 1e-8 * max(1.0, abs(brute.energy)):
                    raise AssertionError(
                        f"Mismatch for x={x}, k={k}, s={s}: alg={alg}, brute={brute}"
                    )
    print(f"All brute-force checks passed for n <= {max_n}, trials={trials}, seed={seed}.")


def parse_points(text: str) -> List[float]:
    return [float(item) for item in text.replace(",", " ").split()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Minimize one-dimensional fixed-cardinality Riesz s-energy by min-cut."
    )
    parser.add_argument("--points", type=str, help="Point coordinates, e.g. '0 1 2 4 7 11'.")
    parser.add_argument("--k", type=int, help="Cardinality of selected subset.")
    parser.add_argument("--s", type=float, default=1.0, help="Riesz exponent; default: 1.")
    parser.add_argument("--examples", action="store_true", help="Run built-in examples.")
    parser.add_argument("--self-test", action="store_true", help="Run randomized brute-force checks.")
    args = parser.parse_args()

    if args.examples:
        run_examples()
    if args.self_test:
        self_test()
    if args.points is not None:
        if args.k is None:
            parser.error("--k is required when --points is given")
        result = minimize_riesz_1d(parse_points(args.points), args.k, args.s)
        print(f"indices_0_based = {result.indices}")
        print(f"indices_1_based = {[i + 1 for i in result.indices]}")
        print(f"points = {result.points}")
        print(f"energy = {result.energy:.12g}")
    if not (args.examples or args.self_test or args.points is not None):
        parser.print_help()


if __name__ == "__main__":
    main()
