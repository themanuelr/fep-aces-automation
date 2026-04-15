# FEP-ACES Automation

This repository automates the setup and analysis of relative binding free energy calculations using the ACES replica-exchange protocol in AMBER. It focuses on ligand transformations on a single protein target, from initial campaign setup to simulation submission and analysis.

## What this repo does

- Organizes a full FEP “campaign” directory (compounds, parametrization, leap, MD, analysis)
- Automates ligand parametrization (ORCA or AM1-BCC/GAFF2)
- Builds all AMBER/ACES input files and job scripts (including λ-window FEP inputs)
- Provides utilities to analyze the resulting FEP trajectories (ΔG/ΔΔG, convergence, RMSD, trajectory conversion)

The code is structured around:

- `FEP_setup.ipynb` — interactive notebook to set up one or more ligand transformations
- `fepfilesetup.py` — core setup logic called by the notebook
- `fepanalysis.py` + `FEP_analysis.ipynb` — analysis utilities and example analysis workflow
- `templates/` — shell scripts and AMBER/ORCA input templates used as blueprints

## Installation and environment

This project assumes:

- A working AMBER/AmberTools installation (for `tleap`, `sander`, `cpptraj`, `antechamber`, `parmchk2`)
- ORCA available 
- Python 3 with `numpy`, `scipy`, `pandas`, `matplotlib` and the usual scientific stack

You will typically:

1. Clone the repository to a location where you can also access your AMBER and ORCA installations.
2. Make sure your environment modules / conda envs for AMBER and ORCA are loaded before running parametrization and MD.

For more detail, see:

- `FEP_setup.ipynb` for the interactive setup sequence.
- `FEP_analysis.ipynb` for example analysis and plotting.

For even more detail, see:

- `docs/FEP_setup_workflow.md` for the full setup workflow.
- `docs/FEP_analysis_workflow.md` for the full analysis workflow.

