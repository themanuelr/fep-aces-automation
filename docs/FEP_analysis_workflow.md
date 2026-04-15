***

## What `.out` Files Contain

Every AMBER MD run produces a `.out` file — a plain-text log written every `NTPR` steps (set to 5,000 in the notebook). Critically for this analysis, each entry contains a line like :

```
 DV/DL  =    Total dU/dl:  0.25   ...   -12.4573
```

This is **∂V/∂λ** — the derivative of the potential energy with respect to the coupling parameter λ at that snapshot. It is the fundamental quantity of **thermodynamic integration (TI)**. Each of the 25 λ-windows produces its own `.out` file, so across all windows you get the full ∂V/∂λ curve from λ=0 (old compound, fully coupled) to λ=1 (new compound, fully coupled).

***

## `process_lines()` — The Core Parser

This is the only function that directly reads the `.out` files.  For a single file (one λ-window), it:

1. Scans every line looking for `'Total dU/dl'`
2. Splits the line and extracts the numeric value at position `[6]` (with a fallback parser for formatting edge cases)
3. Reads the λ value from position `[2]` on the first matching line
4. Returns: `Lambda`, `mean(dU/dl)`, `std(dU/dl)`, and the full raw list `dU_dl[]`
```
"Total dU/dl:  0.25   ...   -12.4573"
                ↑λ               ↑∂V/∂λ value appended each NTPR steps
```


***

## `get_data()` — Assembles All 25 Windows

Calls `get_files()` to glob all `*md.out` files in `03-md/{comp}/{fep_type}/{rep}/`, sorts them by window number, then calls `process_lines()` on each.  Returns four parallel arrays:


| Array | Contents |
| :-- | :-- |
| `LAMBDAS` | The λ value for each window (e.g., 0.0, 0.1, 0.15…1.0) |
| `MEAN_dU_dl` | Mean ∂V/∂λ per window across all simulation time |
| `STD_dU_dl` | Standard deviation per window (spread/noise measure) |
| `ALL_dU_dl` | Full time-series of every single ∂V/∂λ data point per window |


***

## `calc_dG()` — The Central Thermodynamic Integration

Takes `LAMBDAS` and `MEAN_dU_dl` and computes the **free energy change** using the trapezoidal rule:

$$
\Delta G = \int_0^1 \left\langle \frac{\partial V}{\partial \lambda} \right\rangle d\lambda \approx \sum_i \frac{\langle\partial V/\partial\lambda\rangle_i + \langle\partial V/\partial\lambda\rangle_{i+1}}{2} \cdot (\lambda_{i+1} - \lambda_i)
$$

This is the area under the ∂V/∂λ curve. It uses `scipy.integrate.trapezoid()` directly on the mean values. The result is a single **ΔG in kcal/mol** for that compound-replica-fep_type combination.

***

## `get_ddG_table()` — The Final Scientific Output

This is the end goal of the entire pipeline.  For every compound and every replica it calls `calc_dG()`, then:

1. Builds a pandas DataFrame of all ΔG values (compounds as columns, replicas as rows)
2. **Subtracts the unbound ΔG from each replica**: `df[col] = df[col] - df.loc['unbound', col]`
— this converts raw ΔG values into **ΔΔG = ΔG_bound − ΔG_unbound**, which is the **relative binding free energy** between the old and new compound
3. Drops the `unbound` row (already consumed)
4. Computes the mean and standard deviation across replicas as the final reported value, formatted as e.g. `+2.3 (0.4)` (mean ± SD)
5. Prints and returns the formatted DataFrame

The ΔΔG is the number the drug discovery campaign actually cares about: a negative value means the new compound binds **more tightly** to the protein than the old one.

***

## Validation \& Diagnostic Functions

Beyond the core ΔG calculation, the script uses the raw `ALL_dU_dl` time-series data for quality checks:

- **`plot_timeseries()`** — plots the time evolution of ∂V/∂λ for every λ-window as stacked panels; used to spot poor equilibration or instability in specific windows. A red dashed line can be drawn at a `cutoff` time point.
- **`plot_convergence()`** — calculates ΔG on progressively longer fractions of the trajectory (first 10%, 20%…100%) and plots how ΔG evolves; a truly converged simulation produces a flat line.
- **`remove_cutoff(data, cutoff)`** — strips an initial fraction of each window's data (e.g., 5%) before computing means, discarding the equilibration period from the ΔG integral.
- **`calc_running_mean()`** — smooths the raw ∂V/∂λ time series with a sliding window convolution for visual clarity in plots.
- **RMSD functions** — use `cpptraj` to compute backbone RMSD of the protein and ligand RMSD (no-fit, relative to first frame and relative to previous frame), stored in `04-analysis/{comp}/{comp}_{rep}_rmsd/`, to verify the protein and ligand stayed structurally stable throughout the FEP run.
- **Trajectory conversion** — `convert_traj()` writes and executes `cpptraj` scripts to convert AMBER `.crd` binary trajectories into `.xtc` format (GROMACS-compatible), applying `autoimage` to re-center the system across periodic boundaries, stored in `04-analysis/`.

