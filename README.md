# fep-aces-automation
This repository contains a Python-based framework developed during my Postdoctoral work to automate Free Energy Perturbation (FEP) workflows.

Key Features:

    fepfilesetup: Automates the generation of topology and coordinate files, mapping "new" ligands onto an equilibrated "old" reference structure (structural synchronization).

    fepanalysis: A data-processing pipeline that extracts dU/dλ from Amber/ACES output logs, performs numerical integration (TI), and calculates statistical uncertainty.

    HPC Integration: Includes templates and logic for orchestrating 24-window replica exchange simulations on Slurm-managed clusters.

Repository Structure:

    fepfilesetup.py: Core logic for file preparation and template substitution.

    fepanalysis.py: Tools for thermodynamic integration and RMSD validation.

    templates/: Bash and input templates used for automated job submission.

Note: This is research code developed for a specific project. It is shared here to demonstrate the methodology and automation logic used in my published work.
