# Riesz \(s\)-Energy Subset Selection on Ordered Point Sets

This repository contains reproducibility material for the paper

> Michael T. M. Emmerich, *Polynomial-Time Riesz-Energy Subset Selection for Ordered Point Sets on Lines and \(\ell_1\)-Staircases*, arXiv preprint, 2026.

The central problem is the fixed-cardinality minimum Riesz \(s\)-energy subset-selection problem on an ordered point set
\[
x_1 < x_2 < \cdots < x_n.
\]
For a fixed exponent \(s>0\) and cardinality \(k\), the task is to choose
\[
1 \le i_1 < \cdots < i_k \le n
\]
minimizing
\[
E_s(i_1,\ldots,i_k)
  = \sum_{1 \le p < q \le k} (x_{i_q}-x_{i_p})^{-s}.
\]

## Two Python main files

The repository is organized around two executable Python entry points.

```text
riesz_min_cut_main.py              # exact threshold/min-cut solver and small demo
riesz_break_even_reproduction.py   # runtime experiment and enumeration comparison
```

### 1. `riesz_min_cut_main.py`

This is the compact main implementation of the exact one-dimensional threshold/min-cut construction. It provides importable functions for

- constructing deterministic ordered test instances,
- computing Riesz \(s\)-energy values,
- solving the fixed-cardinality problem by the min-cut graph construction,
- checking small instances by complete enumeration.

Run

```bash
python riesz_min_cut_main.py
```

to execute a small deterministic demo with \(n=16\), \(k=8\), \(s=1\). The script prints the selected subset, its energy, the threshold-vector representation, the graph size, and an enumeration check.

### 2. `riesz_break_even_reproduction.py`

This is the empirical runtime reproducibility script used for the manuscript subsection **Empirical Runtime Efficiency**. It compares

1. complete enumeration over all \(k\)-subsets, and
2. the explicit threshold/min-cut construction,

on balanced instances \(n=2k\).

Run

```bash
python riesz_break_even_reproduction.py
```

to create

```text
riesz_break_even_benchmark.csv
riesz_break_even_times.png
```

The benchmark uses deterministic mildly irregular ordered point sets

```python
x_i = i + 0.05 * sin(1.7*i) + 0.01 * (i/n)**2
```

with exponent

```text
s = 1
```

and balanced cardinalities

```text
n = 2k.
```

Complete enumeration is run only up to the range where it remains practical. The min-cut method is run beyond this range.

## Python environment

The scripts were prepared for Python 3.10 or newer. Install the required packages with

```bash
python -m pip install numpy pandas matplotlib networkx numba
```

The algorithm demo requires only `numpy` and `networkx`. The benchmark additionally uses `numba`, `pandas`, and `matplotlib`.

The complete-enumeration baseline in the benchmark is intentionally quite favorable to brute force because it is compiled with Numba. Conversely, the min-cut implementation uses a general-purpose NetworkX backend rather than a specialized max-flow implementation. The reported timings should therefore be interpreted as a reproducibility check and practical break-even experiment, not as an optimized performance claim.

## What is compared in the benchmark?

### Complete enumeration

The exhaustive baseline evaluates all

\[
\binom{n}{k}
\]

candidate subsets and computes their Riesz energy. Since a direct energy evaluation costs \(O(k^2)\), the natural work scale is

\[
k^2 \binom{n}{k}.
\]

### Threshold min-cut construction

The min-cut algorithm uses the threshold-variable graph described in the paper. For \(m=n-k\), the graph has

\[
N = k(n-k) = km
\]

threshold nodes and, in the direct construction, up to

\[
O(k^2(n-k)^2)
\]

finite pairwise arcs. The conservative cut-step bound quoted in the paper is

\[
O(k^4(n-k)^4),
\]

which becomes \(O(k^8)\), equivalently \(O(n^8/256)\), on balanced instances \(n=2k\). A looser comparison scale sometimes written as \(n^8/2\) gives a later theoretical crossing.

## Expected benchmark results

On the benchmark run used for the manuscript, exhaustive enumeration and min-cut agreed on all instances where both were executed. The practical break-even was observed around

```text
n = 24--26,  k = 12--13.
```

A representative timing table is:

| n | k | binomial(n,k) | enumeration time (s) | min-cut time (s) | enumeration / min-cut |
|---:|---:|---:|---:|---:|---:|
| 20 | 10 | 184756 | 0.0089 | 0.0233 | 0.38 |
| 22 | 11 | 705432 | 0.0435 | 0.0404 | 1.08 |
| 24 | 12 | 2704156 | 0.1756 | 0.1827 | 0.96 |
| 26 | 13 | 10400600 | 0.8130 | 0.0787 | 10.33 |
| 28 | 14 | 40116600 | 3.4782 | 0.1106 | 31.44 |
| 30 | 15 | 155117520 | 15.2084 | 0.1718 | 88.53 |

The exact timings depend on the machine, Python version, BLAS/runtime environment, and NetworkX/Numba versions. The main qualitative observation is that the exponential growth of complete enumeration overtakes the polynomial min-cut approach already at moderate balanced sizes.

## Compute environment used for the reported benchmark

The timing values reported above were obtained in a ChatGPT/Python sandbox environment. The precise numbers should be treated as machine-dependent, but the environment was approximately:

```text
Python: 3.13.5
Platform: Linux-4.4.0-x86_64-with-glibc2.41
Processor: unknown
Memory: MemTotal:        4194304 kB
NumPy: 2.3.5
Pandas: 2.2.3
Matplotlib: 3.10.8
NetworkX: 3.6.1
Numba: 0.65.1
```

For publication-quality performance claims, rerun the script on the target machine and report the local environment together with the generated CSV file.

## Theoretical break-even scales

For balanced instances \(n=2k\), comparing

\[
k^2 \binom{2k}{k}
\]

with the conservative balanced min-cut bound

\[
k^8
\]

gives a theoretical crossing at approximately

```text
k = 13, n = 26.
```

If only the number of subsets \(\binom{2k}{k}\) is compared with \(k^8\), the crossing is later, around

```text
k = 19, n = 38.
```

If one compares against the looser scale \(n^8/2\), the corresponding crossings shift to approximately

```text
n = 36
```

with the \(k^2\) energy-evaluation factor, or

```text
n = 48
```

when counting subsets only.

## Reproducibility notes

- The two Python main files are deterministic.
- The point sets are generated internally by the scripts.
- No external data files are required for the runtime benchmark.
- The benchmark CSV stores the selected subset, objective value, graph size, and CPU times.
- The benchmark plot is generated from the same data written to the CSV file.
- The benchmark script checks that complete enumeration and min-cut give the same objective value whenever both methods are run.
- For fair comparisons, report whether the exhaustive baseline is compiled with Numba and which max-flow backend is used for the min-cut computation.

## Citation

If you use this code or benchmark, please cite the accompanying paper.

```bibtex
@misc{Emmerich2026RieszLine,
  author = {Emmerich, Michael T. M.},
  title = {Polynomial-Time Riesz-Energy Subset Selection for Ordered Point Sets on Lines and \(\ell_1\)-Staircases},
  year = {2026},
  eprint = {2606.16946},
  archivePrefix = {arXiv},
  primaryClass = {cs.CG}
}
```

## License

Under CC 4.0 license. All rights reserved.
