# One-Dimensional Minimum Riesz Energy by Min-Cut

This repository contains a small, dependency-free Python implementation of a polynomial-time algorithm for the one-dimensional fixed-cardinality minimum Riesz `s`-energy subset problem.

Given distinct real points

```text
x_0 < x_1 < ... < x_{n-1},
```

a cardinality `k`, and an exponent `s > 0`, the algorithm chooses indices

```text
i_0 < i_1 < ... < i_{k-1}
```

that minimize

```text
sum_{p<q} (x[i_q] - x[i_p])^(-s).
```

The implementation uses the sorted-index lattice and threshold-variable min-cut formulation. It builds one directed `s`-`t` graph and solves one minimum cut.

## Files

```text
riesz_1d_min_cut.py   implementation, examples, command-line interface, self-test
README.md             this file
```

## Requirements

Python 3.10 or newer. No third-party packages are required.

## Quick start

Run the built-in examples:

```bash
python riesz_1d_min_cut.py --examples
```

Run randomized correctness checks against brute force for small instances:

```bash
python riesz_1d_min_cut.py --self-test
```

Solve your own instance:

```bash
python riesz_1d_min_cut.py --points "0 1 2 4 7 11" --k 4 --s 1.5
```

Expected output for this instance:

```text
indices_0_based = [0, 3, 4, 5]
indices_1_based = [1, 4, 5, 6]
points = [0.0, 4.0, 7.0, 11.0]
energy = 0.577850061395
```

## Use as a module

```python
from riesz_1d_min_cut import minimize_riesz_1d

x = [0, 1, 2, 4, 7, 11]
result = minimize_riesz_1d(x, k=4, s=1.5)

print(result.indices)  # zero-based indices in the sorted input
print(result.points)   # selected coordinates
print(result.energy)   # Riesz s-energy
```

## Mathematical idea

Write a feasible `k`-subset as an increasing index vector

```text
i_0 < i_1 < ... < i_{k-1}.
```

Set

```text
i_r = r + y_r.
```

Then feasibility is equivalent to

```text
0 <= y_0 <= y_1 <= ... <= y_{k-1} <= n-k.
```

Introduce threshold variables

```text
z_{r,t} = 1[y_r >= t],     t = 1, ..., n-k.
```

The monotonicity constraints become closure constraints:

```text
z_{r,t+1} <= z_{r,t},
z_{r,t}   <= z_{r+1,t}.
```

For each pair of ranks `p < q`, the pair interaction

```text
V_{pq}(a,b) = (x[q+b] - x[p+a])^(-s)
```

has nonpositive mixed second differences on the feasible threshold domain. This is the one-dimensional Monge inequality for the decreasing convex function `h(d)=d^(-s)`. Hence the quadratic threshold terms are graph-representable:

```text
-w x y = -w x + w x(1-y),     w >= 0.
```

The term `w x(1-y)` is represented by a directed edge `x -> y` of capacity `w`. Unary terms are represented by source/sink arcs. Therefore the Riesz energy equals the cut capacity up to an additive constant, and a minimum cut gives a global optimum.

## Complexity

Let

```text
m = n-k.
```

The graph has

```text
N = k(n-k)
```

nonterminal threshold nodes and

```text
M = O(k^2(n-k)^2)
```

finite pairwise arcs, plus `O(k(n-k))` closure and unary arcs. The implementation solves one max-flow/min-cut problem on this graph.

With a conservative Dinic-type worst-case bound `O(N^2 M)`, the cut step is bounded by

```text
O(k^4 (n-k)^4).
```

For fixed `k`, this gives `O(n^4)` under the same conservative bound. For `k = Theta(n)`, the graph has `O(n^2)` nodes and `O(n^4)` arcs.

## Numerical note

This reference implementation uses floating point capacities. For rational input points and fixed positive integer `s`, all capacities can instead be represented exactly as rational numbers; the same graph construction combined with an exact polynomial-time max-flow algorithm gives a standard Turing-polynomial algorithm.

## Acknowledgement

The sorted-index lattice viewpoint was suggested by Sanchayan Dutta, University of California, Davis, in a MathOverflow discussion. The key point is to view fixed-cardinality feasible sets as a distributive lattice of increasing index vectors, rather than as a cardinality slice of the Boolean lattice.
