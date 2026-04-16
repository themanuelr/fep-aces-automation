## Full Workflow: Step by Step

### Step 1 — `ffs.campaign_dir_setup(campaign_name)`

**Function logic** : Creates the entire campaign directory tree:

```
FEP_campaign/
  00-compounds/rep_struc/
  01-mol_param/
  02-leap/
  03-md/
  04-analysis/
```

Uses `os.makedirs()` for each. Prints instructions to copy 3 starting PDB structures (naming convention: `OLDCOMP_repN_fNNNN.pdb`) into `00-compounds/rep_struc/`.

***

### Step 2 — `ffs.extract_oldcomp_pdb(oldcomp, repN='rep1')`

**Function logic** : Runs a shell `grep` on the starting-structure PDB to extract only lines containing the old compound's residue name, writing `00-compounds/{old}.pdb`. This is the clean ligand-only PDB used as the structural reference.

***

### Step 3 — Manual: Maestro + VMD

No fepfilesetup call. The user opens `{old}.pdb` in **Maestro** (Schrödinger), uses the builder to draw the new compound, saves as `pre_{new}.mol2`, then opens in **VMD** and saves as `{new}.pdb` into `00-compounds/`.

***

### Step 3b — `ffs.get_newcomp_list(oldcomp)`

**Function logic** : `glob`s all `*.mol2` files in `00-compounds/`, strips filenames to 3-letter compound codes, prints and returns `newcomp_list`.

***

### Step 4 — `ffs.generate_cutpdb(oldcomp, newcomp)` [loop over newcomp_list]

**Function logic** :

1. Calls `highlight_newatoms(oldcomp, newcomp, PWD)` — reads both PDB files with `numpy.genfromtxt`, compares atom types **and** coordinates; any atom that differs in type or position is marked as a **softcore (SC) atom**. Returns `new_atoms_SC` and `ori_atoms_SC`.
2. For each replica (`rep1, rep2, rep3`): reads the full starting-structure PDB, writes `{new}_cut.repN.pdb` keeping only the common-core atoms (excluding `ori_atoms_SC`).
3. Appends the SC atom list for both compounds to `00-compounds/SC_atoms.log` — this file drives all downstream FEP mask definitions.

**Output files**: `{new}_cut.rep1/2/3.pdb`, `SC_atoms.log`

***

### Step 5a — `ffs.setup_param('orca', newcomp_list, campaign_name)`

**Function logic** : Routes to `setup_param_orca()`:

- `cp templates/param_orca_temp.sh → 01-mol_param/param_orca_{new}.sh`
- `sed -i 's/NEWCOMP/{new}/g'`

**Template used**: [`param_orca_temp.sh`](https://github.com/themanuelr/fep-aces-automation/blob/main/templates/param_orca_temp.sh) — this script calls a chain of Python helpers: `orca_opt_fileprep.py` (geometry optimization input), `orca_esp_fileprep.py` (ESP input), `orca_resp_fileprep.py` (RESP fitting input), `orca_antechamber_fileprep.py` and `orca_antechamber_badchrg.py` (GAFF typing via antechamber).
User must run `./param_orca_{new}.sh` manually on an HPC (takes ~8 hours).

**Note**: At the end of this file you will find a more detailed description of the ORCA parametrization step.

***

### Step 5b — `ffs.setup_param('am1bcc', ...)` [alternative]

**Function logic** : Routes to `setup_param_am1bcc()`:

- `makedirs 01-mol_param/{new}/`
- `cp templates/am1bcc_temp.sh → am1bcc_{new}.sh`
- `sed NEWCOMP` and `sed CHARGE`
- If `execute_param=True`, runs it in-place via `os.system`

**Template used**: [`am1bcc_temp.sh`](https://github.com/themanuelr/fep-aces-automation/blob/main/templates/am1bcc_temp.sh) — calls `antechamber` (AM1-BCC charge method) and `parmchk2` to generate `{new}.mol2` and `{new}.frcmod`.

***

### Step 6a — `ffs.tleap_setup_organize_dirs(newcomp_list, oldcomp, campaign_name)`

**Function logic** : For every `newcomp × replica (rep1, rep2, rep3, unbound)`:

1. Calls `get_oldcomp_endstring(oldcomp, replica)` — reads starting PDB to find the `TER` line for the old compound (used to know where to insert newcomp in the PDB).
2. `makedirs 02-leap/{new}/{rep}/`
3. `cp templates/get_files_temp.sh → get_files.sh`
4. `sed -i` replaces 7 placeholders: `OLD_MOL`, `STARTFRAMEPATH`, `STARTFRAMENAME`, `NEWMOL`, `REPN`, `ENDSTRING`, `TEMPDIR`
5. Executes `get_files.sh` immediately.

**Template used**: [`get_files_temp.sh`](https://github.com/themanuelr/fep-aces-automation/blob/main/templates/get_files_temp.sh) — this script assembles the tleap working directory: copies `.frcmod` and `.lib` from parametrization, copies and inserts the newcomp coordinates into the protein PDB at the correct `TER` line, and builds `tleap.in` from either `temp_tleap.in` (bound replicas) or `temp_tleap.unbound.in` (unbound replica), plus `temp_getpdb.cpptraj`.

***

### Step 6b — `ffs.tleap_execute(newcomp_list, oldcomp, campaign_name)`

**Function logic** :

- For `unbound`: greps old-compound and newcomp lines from the `rep1` PDB to build an unbound-state PDB.
- Runs `tleap -f tleap.in > tleap.log` → produces `{name}.top` and `{name}.rst`
- Runs `cpptraj getpdb.cpptraj > getpdb.log` → produces `{name}.pdb`

**Output per `02-leap/{new}/{rep}/`**: `.top`, `.rst`, `.pdb` (protein+ligand topology and coordinates for each replica)

***

### Step 7a — `ffs.mdfep_organize_dirs(newcomp_list, campaign_name)`

**Function logic** :

- Creates `03-md/{new}/stdr/` and `03-md/{new}/aces/` with `rep1, rep2, rep3, unbound` subdirectories
- For `unbound/stdr`: also creates `mdsolv/` subdirectory (solvation pre-equilibration)
- Copies `.top` and `.rst` from `02-leap` to the correct `stdr/{rep}/` directories
- Copies `.top` only to `aces/{rep}/` directories
- Renames bound `.rst` files to `00-hNKCC1_{new}.md.rst` (naming convention for simulation restarts)

***

### Step 7b — `ffs.prepare_FEP_files(newcomp, oldcomp, campaign_name, lambda_gap, PROPS_stdr, PROPS_acesR, PROPS_acesUB)` [loop]

**Function logic**  — the most complex step, calls 5 internal functions:

1. **`get_FEP_masks(comp, old_comp)`**: Reads `SC_atoms.log`, parses the softcore atoms for both compounds, returns `FEP_masks` dict with `ini_rname`, `tgt_rname`, `ini_scmask`, `tgt_scmask`.
2. **`produce_input_file(fep_type='stdr', ...)` for rep1/rep2/rep3/unbound**:
    - `extract_clambdas(lambda_gap)` — if `lambda_gap` is a float, builds an evenly spaced list; if a list, uses it directly (e.g., 25 custom λ windows: 0.0, 0.1, 0.150…1.0)
    - `extract_props(PROPS_stdr, 'stdr')` — unpacks `NTPR, NSTLIM, DT, NTC, GTI_BAT_SC`
    - For each λ-window `NN`: `cp templates/temp-FEPstdr-md.in → NN-md.in` and `cp templates/temp-FEPstdr-min.in → NN-min.in`
    - `sed -i` substitutes 10 placeholders: `NTPR, NSTLIM, DT, NTC, GTI_BAT_SC, CLAMBDA, TIMASK1, TIMASK2, SCMASK1, SCMASK2, N_WINDOWS, ALL_WINDOWS`
    - **Unbound solvation special case**: also copies `temp-FEPstdr-01-min.in`, `temp-FEPstdr-02-md.in`, `temp-FEPstdr-03-md.in`, `temp-FEPstdr-04-md.in` into `mdsolv/` and applies the same seds.
3. **`produce_input_file(fep_type='aces', ...)` for rep1/rep2/rep3**:
    - `extract_props(PROPS_acesR, 'aces')` — also unpacks `NUMEXCHG`
    - `cp templates/temp-FEPaces-md.in → NN-md.in` per λ-window
    - Same sed substitutions + `NUMEXCHG`
4. **`produce_input_file(fep_type='aces', ...)` for unbound**:
    - Same as above but uses `PROPS_acesUB` (3× longer `NSTLIM`)
5. **`create_groupfile(comp, 'aces', ...)` × 2** (once for `[rep1,rep2,rep3]`, once for `[unbound]`):
    - Globs all `*md.in` files in `aces/rep1/` to get window numbers
    - Writes `groupfile_windowpart_0` listing every window's AMBER command line: `-O -p {TOP} -i NN-md.in -c ../../stdr/{rep}/NN-{name}.FEPstdr.md.rst -r … -x … -o … -inf …`
    - This groupfile drives the ACES replica-exchange protocol across all 25 λ-windows simultaneously.

***

### Step 7c — `ffs.distribute_submission_jobs(newcomp_list, campaign_name)`

**Function logic** : For both `stdr` and `aces`, for every compound and every replica:

- `cp templates/temp-job_FEPstdr_1.sh → job_FEPstdr_1.sh`
- `cp templates/temp-job_FEPaces_1.sh → job_FEPaces_1.sh`
- `sed` replaces `COMP, REPN, NREP, WIN, GROUPFILE, PREFIX`
- For unbound `stdr`: also `cp templates/temp-job_mdsolv_1.sh → mdsolv/job_mdsolv_1.sh`

**HPC Submission Chain** :

```
qsub job_mdsolv_1.sh          ← user submits only this one
  └─ submits all stdr job_FEPstdr_1.sh (rep1,2,3,unbound)
       └─ each stdr job submits its aces job_FEPaces_1.sh
            └─ ACES runs replica exchange across 25 λ-windows
```


***

## Complete Template Inventory

| Template File | Used By | Purpose |
| :-- | :-- | :-- |
| `param_orca_temp.sh` | `setup_param_orca()` | ORCA QM parametrization script |
| `am1bcc_temp.sh` | `setup_param_am1bcc()` | AM1-BCC charge + GAFF parametrization |
| `orca_opt_fileprep.py` | `param_orca_temp.sh` | Geometry optimization input |
| `orca_esp_fileprep.py` | `param_orca_temp.sh` | ESP calculation input |
| `orca_resp_fileprep.py` | `param_orca_temp.sh` | RESP charge fitting input |
| `orca_antechamber_fileprep.py` | `param_orca_temp.sh` | Antechamber GAFF typing |
| `orca_antechamber_badchrg.py` | `param_orca_temp.sh` | Charge correction helper |
| `get_files_temp.sh` | `tleap_setup_organize_dirs()` | Assembles tleap working dir |
| `temp_tleap.in` | `get_files_temp.sh` | tleap input (bound replicas) |
| `temp_tleap.unbound.in` | `get_files_temp.sh` | tleap input (unbound) |
| `temp_getpdb.cpptraj` | `get_files_temp.sh` | cpptraj PDB extraction script |
| `temp-FEPstdr-md.in` | `produce_input_file()` | Standard AMBER FEP MD input |
| `temp-FEPstdr-min.in` | `produce_input_file()` | Standard AMBER FEP minimization |
| `temp-FEPaces-md.in` | `produce_input_file()` | ACES replica-exchange MD input |
| `temp-FEPaces-min.in` | `produce_input_file()` | ACES minimization input |
| `temp-FEPstdr-01-min.in` | `produce_input_file()` (unbound) | Unbound solvation: minimization |
| `temp-FEPstdr-02-md.in` | `produce_input_file()` (unbound) | Unbound solvation: heating |
| `temp-FEPstdr-03-md.in` | `produce_input_file()` (unbound) | Unbound solvation: equilibration |
| `temp-FEPstdr-04-md.in` | `produce_input_file()` (unbound) | Unbound solvation: production prep |
| `temp-job_FEPstdr_1.sh` | `distribute_submission_jobs()` | HPC stdr job script |
| `temp-job_FEPaces_1.sh` | `distribute_submission_jobs()` | HPC aces job script |
| `temp-job_mdsolv_1.sh` | `distribute_submission_jobs()` | HPC mdsolv job script (entry point) |


## Step 5a — in detail: `ffs.setup_param('orca', newcomp_list, campaign_name)` 

## Stage 1: `01-opt/` — Geometry Optimization

**Script**: `orca_opt_fileprep.py`

This script reads `00-compounds/pre_{new}.mol2` — the raw molecule drawn in Maestro. It parses the `@<TRIPOS>ATOM` block to extract every atom's element and (x,y,z) coordinates, then writes an ORCA input file:

```
! B3LYP 6-31g* opt

*xyz 0 1
 C    1.234000    0.567000   -0.890000
 N    2.345000   ...
 ...
*
```

The keywords `B3LYP 6-31g* opt` instruct ORCA to run a **DFT geometry optimization** at the B3LYP/6-31G* level of theory. This finds the true energy-minimum geometry of the molecule in vacuum — the "relaxed" structure that the molecule naturally prefers.

The shell script then runs: `orca opt_{new}.inp > opt_{new}.out`

**Output**: `opt_{new}.xyz` — optimized Cartesian coordinates. This is the structure used for all downstream calculations.

***

## Stage 2: `02-esp/` — ESP Calculation

**Script**: `orca_esp_fileprep.py`

This reads the optimized geometry from `01-opt/opt_{new}.xyz` and writes a new ORCA input:

```
! HF 6-31g* chelpg

*xyz 0 1
 C    1.234000    ...
 ...
*
```

The key here is `HF 6-31g* chelpg`:

- **HF** (Hartree-Fock): the quantum chemistry method used — deliberately chosen at this level because the AMBER/GAFF force field was parameterized against HF/6-31G* charges, so consistency is required.
- **6-31G***: a basis set (the mathematical representation of electron orbitals).
- **chelpg**: instructs ORCA to compute the **CHELPG electrostatic potential** on a grid of points surrounding the molecule (using a Connolly surface), outputting a `.vpot` file with the potential values at every grid point.

The shell script then runs: `orca esp_{new}.inp > esp_{new}.out`

It also grabs the total molecular charge from the output:

```bash
charge=$(grep "Sum of atomic charges:" ../02-esp/esp_{new}.out | awk '{print $5}')
charge_int=$(echo "$charge/1" | bc)
```

This integer charge is passed to every subsequent step.

**Output**: `esp_{new}.scfp.esp_{new}.vpot` — the raw electrostatic potential grid file.

***

## Stage 3: `03-resp/` — RESP Charge Fitting

**Script**: `orca_resp_fileprep.py`

This is a format-conversion step. The `.vpot` file from ORCA is not in the format expected by `pyresp`/`resp`. This script:

1. Reads the `.vpot` file and `opt_{new}.xyz` atom list.
2. Uses `PyAstronomy.pyasl.AtomicNo()` to convert element symbols (e.g., `C`, `N`, `O`) into atomic numbers (6, 7, 8).
3. Rewrites the file in the exact column format required by the RESP tools, outputting `esp_{new}.vpot`.

The shell script then runs the **two-stage RESP fit** :

```bash
pyresp_gen.py -i esp_{new}.vpot -f1 resp1.in -f2 resp2.in -p chg -q $charge_int
py_resp.py -O -i resp1.in -o resp1.out -t q1 -e esp_{new}.vpot -s resp1.esp
py_resp.py -O -i resp2.in -o resp2.out -t q2 -q q1 -e esp_{new}.vpot -s resp2.esp
```

- **Stage 1** (`resp1`): fits all charges with tight restraints — constrains chemically equivalent atoms (e.g., all H atoms on a methyl group get the same charge).
- **Stage 2** (`resp2`): refines buried/aliphatic atom charges with relaxed restraints, using Stage 1 charges as a starting point (`-q q1`).

**Output**: `resp2.out` — contains the final `q(opt)` column with the RESP charge for every atom.

***

## Stage 4: `04-antechamber/` — GAFF Typing + Final mol2

This stage has two scripts working together.

### `orca_antechamber_fileprep.py`

This script solves a coordinate mismatch: the original `{new}.pdb` has the user-drawn Maestro coordinates, but the RESP charges were derived from the **ORCA-optimized geometry**. For the charges to make physical sense, the topology must use the same geometry. This script:

1. Reads `00-compounds/{new}.pdb` for the **PDB atom metadata** (atom names, residue info, record format).
2. Reads `01-opt/opt_{new}.xyz` for the **optimized Cartesian coordinates**.
3. Merges them: rewrites every `ATOM` line keeping PDB metadata but substituting in the ORCA-optimized (x,y,z).
4. Writes `opt_{new}.pdb` — a PDB file with correct atom names AND optimized geometry.

The shell script then runs **antechamber** :

```bash
antechamber -i opt_{new}.pdb -fi pdb -o badchrg_{new}.mol2 -fo mol2 \
            -c gas -nc ${charge_int} -rn ${new} -at gaff2
```

This assigns **GAFF2 atom types** (e.g., `c3`, `n2`, `oh`) to every atom — these are the force field type labels that define which bonded parameters (bonds, angles, dihedrals) to use. The `-c gas` flag tells antechamber to compute its own gas-phase charges (these are `badchrg` — "bad charges" — because they will be **replaced** in the next step).

### `orca_antechamber_badchrg.py`

This is the critical charge-replacement step. Antechamber assigned GAFF2 atom types but used its own charge method. This script:

1. Reads the **RESP charges** from `03-resp/resp2.out`, extracting every `q(opt)` value.
2. Reads `badchrg_{new}.mol2` (the antechamber output with wrong charges).
3. Overwrites the charge column in the `@<TRIPOS>ATOM` section with the RESP-derived values, writing the final `{new}.mol2`.

The shell then runs **parmchk2** :

```bash
parmchk2 -i {new}.mol2 -f mol2 -o {new}.frcmod
```

This generates the `.frcmod` file — a **force field modification file** containing any GAFF2 parameters (bonds, angles, dihedrals, van der Waals) that are missing for atom type combinations specific to this molecule.

Finally, `{new}.mol2` and `{new}.frcmod` are copied to the parent `01-mol_param/{new}/` directory, ready to be consumed by `tleap` in Step 6a.

***

## Full Pipeline Summary

```
pre_{new}.mol2 (Maestro)
       │
       ▼
 01-opt/  orca_opt_fileprep.py
          → opt_{new}.inp  [B3LYP/6-31G* geometry optimization]
          → ORCA runs → opt_{new}.xyz  ← relaxed geometry
       │
       ▼
 02-esp/  orca_esp_fileprep.py
          → esp_{new}.inp  [HF/6-31G* CHELPG ESP calculation]
          → ORCA runs → esp_{new}.vpot  ← electrostatic potential grid
                      → charge_int extracted from esp_{new}.out
       │
       ▼
 03-resp/  orca_resp_fileprep.py  ← reformats .vpot for pyresp
           → pyresp_gen.py generates resp1.in, resp2.in
           → py_resp.py stage 1: tight restraints → q1
           → py_resp.py stage 2: relaxed restraints → resp2.out  ← q(opt) per atom
       │
       ▼
 04-antechamber/  orca_antechamber_fileprep.py
                  → merges PDB metadata + ORCA coords → opt_{new}.pdb
                  → antechamber assigns GAFF2 atom types → badchrg_{new}.mol2
                  orca_antechamber_badchrg.py
                  → replaces charges with RESP q(opt) values
                  → fixes ne→n2 atom type
                  → final {new}.mol2  ← GAFF2 types + RESP charges
                  → parmchk2 → {new}.frcmod  ← missing parameters
       │
       ▼
  {new}.mol2 + {new}.frcmod  →  copied to 01-mol_param/{new}/
  (consumed by tleap in Step 6a)
```
