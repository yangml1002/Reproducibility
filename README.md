# Reproducibility Material for AIPE-TC Numerical Experiments

## Manuscript

Quasi-synchronization of stochastic heterogeneous complex networks via aperiodically intermittent pinning event-triggered control

Manuscript ID: CHAOS-D-26-05548

## Repository contents

```text
Reproducibility/
├── README.md
├── code/
├── data/
└── tables/
```

- `code/`: MATLAB simulation, certificate, post-processing, figure-generation, and table-generation scripts.
- `data/`: precomputed Monte Carlo, certificate, sensitivity, and analytical-comparison data supporting the numerical results. All precomputed data files are stored directly in `data/`.
- `tables/`: human-readable LaTeX/numerical summaries corresponding to the reported numerical tables.

The repository does not store generated figure files. The plotting scripts in `code/` can regenerate the numerical figures from the included precomputed data when needed.

## Numerical environment

Verified locally using MATLAB `version` and `ver`:

- MATLAB 23.2.0.2365128 (R2023b).
- Microsoft Windows 10, version 10.0, build 19045.
- Monte Carlo realizations: 100.
- Random-number generator: Mersenne Twister; seed: 1.
- Euler–Maruyama step: `1e-5 s`.
- Baseline sampling period: `5e-5 s`.
- Baseline horizon: `1.5 s`; the persistent-mismatch experiment uses `5 s`.

Each realization uses a common scalar Brownian increment across all ten nodes of the stochastic second-order oscillator network.

## Main reproduction workflow

Run scripts from `code/` in MATLAB. The scripts resolve data and output directories relative to their own location. All scripts read from or write final numerical outputs directly to `data/`; no nested data subdirectories are required.

### Full experiment drivers (not needed for post-processing)

1. `run_revision_final_AD.m`: main A–D comparisons, threshold grid, and persistent-mismatch experiments.
2. `run_reviewer4_remaining.m`: sparse pinning and one-at-a-time sensitivity experiments.
3. `run_sampling_sweep_consistent.m`: sampling sweep using the main simulation engine and Brownian-path convention.
4. `run_R4_sparse_certified_addon.m`: theorem-certified sparse placement comparison.

These drivers recalculate experiments and overwrite their corresponding result files. To retain the supplied reference data, run full drivers only in a separate working copy. No full Monte Carlo driver was rerun during preparation of this repository.

### Post-processing using precomputed data

Run these individually; no single script is claimed to produce all outputs:

| Script | Inputs | Outputs / purpose |
|---|---|---|
| `figure2upperbounds.m` | Fixed analytical comparison parameters | Fig. 3, bound-comparison CSV and parameter record |
| `make_final_AD_figures.m` | Main/bias MAT results, bias CSVs, schedule boundaries | Figs. 2, 4, 6 and additional diagnostic figures |
| `make_Figure5_tradeoff.m` | `final_main_summary.csv` | Three-panel Fig. 5 |
| `make_R4_sparse_figure.m` | `r4_sparse_summary.csv`, `sparse_certified_addon_summary.csv` | Fig. 7 |
| `make_R4_sampling_figure.m` | `sampling_consistent_summary.csv`, `r4_schedule_checks.csv` | Fig. 8 |
| `make_R4_sensitivity_figure.m` | `r4_sensitivity_summary.csv` | Additional sensitivity curves |
| `make_final_AD_tables.m` | Main/bias summaries, certificates, paired CSV | Tables 3–5, paired-difference CSV, highlights |
| `make_R4_remaining_tables.m` | Sparse/add-on/sampling/sensitivity summaries | Tables 6–7 and full sensitivity table |

The simulation engines and certificate/helper scripts are dependencies, not standalone full experiment drivers.

## Manuscript-symbol mapping

- Internal raw-data fields `E_input` and `E_mean` correspond to the actuation-effort metric denoted by $J_u$ in the revised manuscript.
- Dataset identifiers such as `grid_x0.1` use `x` internally for the threshold multiplier $\xi$.
- Variables named `kappa` and `kappa_matrix` inside the simulation engine denote pinning gains $\kappa_i$, not the threshold multiplier.
- `deltaH` in the Halanay comparison script denotes $\delta_{\rm H}$; the parameter record uses `delta_H`.
- Raw certificate values retain computational precision; the current manuscript displays $h_{\rm cert}=6.998\times10^{-5}\,\mathrm{s}$ for the complete-graph reference configuration.
- The full sensitivity table's “Matrix/sampling conditions” column summarizes the displayed matrix and sampling checks and is not a complete independent re-derivation of every hypothesis of Theorem 1.

## Numerical-zero convention

The raw $q_h=0$ tail-window values are retained in the data files at their computed $O(10^{-35})$ scale. In agreement with the revised manuscript, the human-readable summary reports this case as approximately numerical zero rather than with artificial significant digits. Raw confidence intervals are preserved.

The manuscript positive-mismatch fit uses only $q_h=0.02,0.05,0.10$, giving $R^2=0.999268776\approx0.9993$. The $q_h=0.20$ case is retained only for additional diagnostics and is not included in manuscript figures, tables, or the fit. Additional diagnostic figures may display it.

## Monte Carlo design

The same Brownian paths are reused within paired strategy comparisons. Confidence intervals quantify Monte Carlo sampling variability of the discretized experiment and do not represent time-discretization error.

The refresh at $t=0$ is counted once. Active-start refreshes are compulsory; active-end switch-off and the terminal endpoint add no state update. $N_u$ counts global controller updates; $N_s$ counts sensing/trigger-evaluation opportunities.

$J_\varrho$ is the time-averaged squared synchronization error, evaluated with trapezoidal integration. $J_u$ is the time-averaged sum of squared control inputs accumulated from the actual held input.

## Precomputed data

Precomputed data are included so that manuscript figures and tables can be regenerated without rerunning the full 100-path Monte Carlo simulations. Full drivers are provided for recalculation from scratch.

The analytical all-window schedule deficit is $0.013125$, consistent with the finite-horizon endpoint check. The admissible $N_c=0.05$ and $\lambda_c=0.75$ remain unchanged.
