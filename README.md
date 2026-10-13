# Computational verification — Combinatorics on Words and Symbolic Dynamics (Fibonacci word density)

Code accompanying the PhD thesis of Duaa Abdullah (MIPT, 2026), supervisor A. Ya. Belov.
All finite-computation claims of the thesis are re-checked by a single script.

## Run
```bash
pip install -r requirements.txt     # numpy, scipy, networkx
python verify_all.py                # all sections (~ a few minutes)
python verify_all.py A F            # selected sections
```
Output: `PASS`/`FAIL` per claim; exit code 1 if anything fails. `verification_log.txt` is the log of a full run
(Python 3.12).

## Sections ↔ thesis
| Sec. | Thesis | What is checked |
|---|---|---|
| A | 3.4 | exact count `floor((n+1)/φ²)`, discrepancy bounds and their sharpness, balance of all factors, local structure, p(n)=n+1 |
| B | 3.4 | numerical examples/tables (381, 145 occurrences, palindrome tables, s_20, error formula) |
| C | 3.5 | Beatty description of the positions of 0 and 1 |
| D | 3.5 | framed words q_m: exact error formula, table |
| E | 3.5 | Millin-type telescoping series, d'Ocagne, symmetry identity, sum over 1/(F_k F_{k+2r}) |
| F | 3.6 | palindromic prefixes (Fibonacci, ternary y, characteristic Sturmian, m-bonacci), richness of y_n, density values |
| G | 3.8 | generalized Thue–Morse: κ(N,b,j) vs brute force, b=2 corollary, complexity observation |
| H | 5.1 | factor trees: |V|, irr, irr₂, M₁, M₂, Table for T_8 |
| I | 5.1 | all trees n≤14: Prop., conjecture irr₂<irr₂(star); Lemma on adding an edge (random graphs) |
| J | 5.1 | 2-D recursion T(u,v) |
| K | 2.2 | sanity check of the n! word (finite; does not prove the theorem) |
| L | 4.3 | Welch t, CI, d, power table, η² recomputed from the printed summary statistics |

## Stand-alone scripts (`scripts/`)
Small independent scripts (one per computation in the thesis); run e.g. `python scripts/03_factor_tree_indices.py`.

| Script | Thesis |
|---|---|
| `01_factor_complexity_W.py` | Observation on p(k) of W (3.8) |
| `02_repetitions_search.py` | computer search for (3,σ^j)-repetitions, N=3..6, block length < 400 |
| `03_factor_tree_indices.py` | Theorem on factor trees, table of T_8 |
| `04_conjecture_irr2_trees.py` | conjecture irr₂(T) < irr₂(S_n), all trees n ≤ 12 |
| `05_prefix_contains_all_factors.py` | which s_n contain all factors of length ≤ 30; palindromes per length |
| `06_kappa_bruteforce.py` | κ(N,b,j) versus brute force |
| `07_density_101_plot.py` | convergence of the density of 101 |
| `08_factor_frequencies.py` | frequencies of the length-3 factors (101 → 1/φ⁴) |
| `09_figure_letter_density.py` | figure `figfib_conv.pdf` |
| `10_palindromes_in_prefixes.py` | table of palindrome densities in s_4, s_5, s_6 |
| `11_palindromic_subsequences.py` | palindromic factors vs. scattered palindromes |

Convention: everywhere s_n = Θⁿ(0) with Θ: 0→01, 1→0 (so the rare letter is 1).

## Light versions for online compilers (`online_light/`)
Three short scripts (< 1 s each) for online Python sandboxes with a 3-second limit:
`online_1_fibonacci.py`, `online_2_trees_thuemorse.py`, `online_3_kappa.py`.

## Not covered
* ANCOVA, paired t, Mann–Whitney, Shapiro–Wilk of Section 4.3 need the raw pupil data (anonymised CSV).
  Put it in `data/` and add `analysis_ancova.py`.
* Theorems proved on paper (all-n statements) are only tested on finite ranges; a PASS is evidence, not a proof.

## Citation
Please cite the thesis and the archived release (Zenodo DOI: *to be added*).
