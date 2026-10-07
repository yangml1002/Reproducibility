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

Generated figure files are not stored in the repository. The final A--D plotting script and the targeted persistent-mismatch post-processor write to a temporary output folder outside the repository; this disposable output is not part of the archived reproducibility material.

## Numerical environment

Verified locally using MATLAB `version` and `ver`:

- MATLAB 23.2.0.2365128 (R2023b).
- Microsoft Windows 10, version 10.0, build 19045.
- Monte Carlo realizations: 100.
- Random-number generator: Mersenne Twister; seed: 1.
- Euler–Maruyama step: `1e-5 s`.
- Baseline sampling period: `5e-5 s`.
- Baseline horizon: `1.5 s`; the persistent-mismatch experiment uses `5 s`.
- Target system: `d0 = 0.5` and `zeta_star(0) = (0,0)^T`, so the target trajectory remains identically zero. The numerical implementation therefore evolves the node/error states directly relative to this zero target rather than integrating a separate target trajectory.

Each realization uses a common scalar Brownian increment across all ten nodes of the stochastic second-order oscillator network.

## Main reproduction workflow

Run scripts from `code/` in MATLAB. The scripts resolve data and output directories relative to their own location. All scripts read from or write final numerical outputs directly to `data/`; no nested data subdirectories are required.

### Full experiment drivers (not needed for post-processing)

1. `run_revision_final_AD.m`: main A–D comparisons, threshold grid, and persistent-mismatch experiments.
2. `run_reviewer4_remaining.m`: legacy Reviewer-4 sparse/preliminary diagnostic driver; its preliminary sensitivity branch is not the source of the finalized Table 8. Finalized one-at-a-time sensitivity results are produced by `run_parameter_sensitivity.m`. The legacy driver does not recreate `r4_sensitivity_summary.csv`.
3. `run_sampling_sweep_consistent.m`: sampling sweep using the main simulation engine and Brownian-path convention.
4. `run_R4_sparse_certified_addon.m`: theorem-certified sparse placement comparison.

These drivers recalculate experiments and overwrite their corresponding result files. To retain the supplied reference data, run full drivers only in a separate working copy. Results for Tables 8–10 are supplied as precomputed reference data; their dedicated drivers are described below.

### Post-processing using precomputed data

Run these individually; no single script is claimed to produce all outputs:

| Script | Inputs | Outputs / purpose |
|---|---|---|
| `figure2upperbounds.m` | Fixed analytical comparison parameters | Fig. 3, bound-comparison CSV and parameter record |
| `make_final_AD_figures.m` | Main/bias MAT results, bias CSVs, schedule boundaries | Figs. 2, 4, 6 and additional diagnostic figures |
| `recompute_certificates.m` | Existing case metadata in precomputed MAT/CSV files | Current certificate CSVs and `certificate_bound_comparison.csv`; deterministic verification before publication |
| `update_mismatch_certificate_table.m` | Bias summary and `certificate_bound_comparison.csv` | Updates only the certificate column of Table 5; emits a regression report in the specified work folder |
| `make_persistent_mismatch_figure.m` | Existing bias curves, summary, and current certificates | Fig. 6 only; no simulation, with optional output-folder argument |
| `make_Figure5_tradeoff.m` | `final_main_summary.csv` | Manuscript Fig. 5; `threshold_tradeoff_2panel.pdf` (three panels) |
| `make_R4_sparse_figure.m` | `r4_sparse_summary.csv`, `sparse_certified_addon_summary.csv` | Fig. 7 |
| `make_R4_sampling_figure.m` | `sampling_consistent_summary.csv`, `r4_schedule_checks.csv` | Fig. 8 |
| `make_R4_sensitivity_figure.m` | `parameter_sensitivity_full.csv` | Optional finalized sensitivity curves; output goes to a temporary folder outside the repository |
| `make_final_AD_tables.m` | Main/bias summaries, certificates, paired CSV | Tables 3–5, paired-difference CSV, highlights |
| `make_R4_remaining_tables.m` | Sparse/add-on/sampling summaries | Tables 6–7 only |

The simulation engines and certificate/helper scripts are dependencies, not standalone full experiment drivers.

## Manuscript-symbol mapping

- Internal raw-data fields `E_input` and `E_mean` correspond to the actuation-effort metric denoted by $J_u$ in the revised manuscript.
- Dataset identifiers such as `grid_x0.1` use `x` internally for the threshold multiplier $\xi$.
- Variables named `kappa` and `kappa_matrix` inside the simulation engine denote pinning gains $\kappa_i$, not the threshold multiplier.
- `deltaH` in the Halanay comparison script denotes $\delta_{\rm H}$; the parameter record uses `delta_H`.
- Raw certificate values retain computational precision; the current manuscript displays $h_{\rm cert}=6.998\times10^{-5}\,\mathrm{s}$ for the complete-graph reference configuration.
- Table 8 uses the finalized one-at-a-time data in `parameter_sensitivity_full.csv` and is generated by `make_tables_8_10.m`. Its certificate status combines the average-control-rate all-window check with the active/rest matrix and sampling conditions.

## Numerical-zero convention

The raw $q_h=0$ tail-window values are retained at their computed $O(10^{-35})$ scale. In the human-readable manuscript table, this case is reported as approximately zero (numerical), while the raw confidence intervals are preserved.

The manuscript positive-mismatch fit uses only $q_h=0.02,0.05,0.10$, giving $R^2=0.999268776\approx0.9993$. The $q_h=0.20$ case is retained only for additional diagnostics and is not included in manuscript figures, tables, or the fit. Additional diagnostic figures may display it.

## Monte Carlo design

The same Brownian paths are reused within paired strategy comparisons. Confidence intervals quantify Monte Carlo sampling variability of the discretized experiment and do not represent time-discretization error.

The refresh at $t=0$ is counted once. Active-start refreshes are compulsory; active-end switch-off and the terminal endpoint add no state update. $N_u$ counts global controller updates; $N_s$ counts sensing/trigger-evaluation opportunities.

$J_\varrho$ is the time-averaged squared synchronization error, evaluated with trapezoidal integration. $J_u$ is the time-averaged sum of squared control inputs accumulated from the actual held input.

## Precomputed data

Precomputed data are included so that manuscript figures and tables can be regenerated without rerunning the full 100-path Monte Carlo simulations. Full drivers are provided for recalculation from scratch.

The analytical all-window schedule deficit is $0.013125$, consistent with the finite-horizon endpoint check. The admissible $N_c=0.05$ and $\lambda_c=0.75$ remain unchanged.

### Current analytical certificates

The current ultimate certificate is $R_{\rm av}$ from Theorem 1. With $d=\min\{r,\Lambda_2\}$, $\bar r=d+\lambda_c(r-\Lambda_2)_+$, and $T_c=N_c/\lambda_c$,

$$R_{\rm av}=2\Lambda_4\left[\frac{1-e^{-dT_c}}d+\frac{e^{-dT_c}}{\bar r}\right].$$

The field `R_ultimate` equals `R_av`, preserving compatibility with downstream readers. The former minimum-rate quantity $R_{\rm minrate}=2\Lambda_4/d$ is retained only as an analytical comparison. The theorem retains the schedule-resolved kernel $\mathcal K$ in the transient bound, including the faster local disturbance propagation inherited from Lemma 8. The closed-form ultimate certificate $R_{\rm av}$ is obtained after the relaxation $\mathcal K\le\mathcal P$ and the subsequent application of the all-window average-control-rate envelope $\Psi$. Thus, the local $\mathcal K$-kernel refinement and the improvement of the closed-form average-rate certificate play distinct roles.

Run `recompute_certificates(work_dir)` to regenerate the certificate CSVs from existing precomputed case metadata without rerunning Monte Carlo simulations. The routine checks matrix tests, sampling limits, decay roots, disturbance constants, and existing statuses before publishing any CSV. `certificate_bound_comparison.csv` compares the average-rate and minimum-rate bounds; rows outside the certified set are algebraic comparisons, not theorem-certified bounds. CP-PETC rows are not certified by this intermittent theorem. The named `certified` field combines the active/rest matrix tests, sampling condition, and independently verified all-window schedule condition.

All precomputed data remain directly under `data/`. Stochastic MAT files and their per-path metrics are preserved unchanged; their embedded certificate arrays are metadata snapshots rather than the authoritative ultimate certificates. Use the named certificate CSVs (`parameter_sensitivity_certificates.csv`, `performance_match_certificates.csv`, and `sparse_certified_addon_certificates.csv` for the corresponding experiments) for the `R_av` values.

The finite-horizon schedule flags are preserved separately from `schedule_all_window_verified`. For the rounded-duty engine, an all-window certificate is recorded only for exactly realizable quarter-fraction duties, whose complete periods have a multiple-of-four step count; other rounded-duty rows remain algebraic sensitivity comparisons. The sensitivity and matched-comparison rows use their verified all-window schedule statuses. A finite-horizon endpoint check alone is not promoted to an all-time proof.

## Tables 8--10: sensitivity, matched comparison, and time-step verification

The finalized additions use the main simulation engine, 100 paired paths, Twister seed 1, and the same baseline counting and metric conventions. Only the reported numerical outputs are included in this repository.

Run these drivers only if recalculation is needed, in this order:
1. run_parameter_sensitivity.m: the 15 points in Table 8, keeping other baseline settings fixed.
2. run_performance_matched_activation_comparison.m: the CP-PETC reference and actual AIPE-TC activation-rate grid; Table 9.
3. run_timestep_refinement_check.m: reads the selected matched workspace and physical schedule; Table 10.
4. make_tables_8_10.m: post-processes supplied finalized CSVs into Tables 8–10 without simulation, reproducing manuscript-style scientific notation and the Table 9 selection note.

For the near-error-matched comparison, both event-trigger thresholds are fixed at xi=1. CP-PETC remains continuously active, while AIPE-TC is evaluated at the actually simulated activation-rate grid 0.9, 0.95, 0.975, 0.9875, and 0.99375. The selected actually simulated point minimizes the symmetric relative mean-error gap. No interpolation is used. The selected rate is 0.99375, with a 0.1383% residual gap. CP-PETC has no intermittent activation-rate parameter (shown as -- in Table 9).

At the selected near-error-matched point, AIPE-TC reduces active-time usage by approximately 0.6173% and sensing opportunities by approximately 0.5767% relative to CP-PETC, while requiring more controller updates and nearly the same actuation effort.

For variable activation rates, the ceiling rule guarantees that active duration in every complete cycle is at least lambda_c times its period. With P_max=0.07 s, the all-window deficit is bounded by lambda_c(1-lambda_c)P_max. The finalized certificate includes this check, the active/rest matrix inequalities, and sampling inequality. Baseline duty=0.75 is unchanged.

Paired 95% confidence intervals use actual pathwise differences AIPE-TC minus CP-PETC, with mean(d) +/- t_(0.975,99)*std(d)/sqrt(100); marginal confidence intervals are not subtracted.

The time-step check uses a coarse-anchored Brownian bridge. Original coarse increments are generated exactly as in Table 9; an independent bridge splits each into two fine increments. Summing each fine pair recovers the original coarse increment to floating-point precision. Coarse and fine simulations share the physical schedule, initial conditions, thresholds, and sampling period. Table 10 coarse rows agree pathwise with Table 9. Temporal discretization is checked separately from Monte Carlo variability.

| Manuscript output | Finalized files directly in data/ | Driver |
|---|---|---|
| Table 8 | parameter_sensitivity_full.csv; parameter_sensitivity_workspace.mat | run_parameter_sensitivity.m |
| Table 9 search | performance_match_activation_grid.csv | run_performance_matched_activation_comparison.m |
| Table 9 statistics | performance_match_summary.csv; performance_match_paired.csv; performance_match_workspace.mat; same_threshold_paired.csv | matching driver; paired_ci.m |
| Table 9 certificate | performance_match_certificates.csv; schedule_all_window_checks.csv; performance_match_certificate.mat (preserved metadata snapshot) | recompute_certificates.m; verify_reproducibility.m |
| Table 10 | timestep_refinement_matched_summary.csv; timestep_refinement_matched_paired.csv; timestep_refinement_matched_workspace.mat | run_timestep_refinement_check.m |
| Baseline identity | baseline_reproduction.mat | finalized drivers |

Workspaces preserve per-path metrics at full precision. Existing main, sparse, mismatch, and sampling datasets remain unchanged. Human-readable tables round displayed statistics only.

`certificate_R4_remaining.m` is a certificate helper, not a standalone driver for generating `schedule_all_window_checks.csv`.

## Verification utilities

`verify_reproducibility.m` uses only the supplied precomputed data and performs no Monte Carlo simulation. It checks:

- sensitivity baseline pathwise identity;
- all-window average-control-rate bounds;
- Table 9 / Table 10 coarse pathwise identity;
- Brownian-bridge aggregation;
- paired confidence intervals;
- near-error-matched gap and resource reductions.

It produces `data/reproducibility_checks.csv` and reconstructs `data/schedule_all_window_checks.csv` after checking consistency with the supplied values.

`smoke_test.m` is an optional two-path, very-short-horizon engine check. It is not part of the manuscript numerical results and does not overwrite the supplied reference data.
