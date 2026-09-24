# Radiation Shielding Research

## Public reproducibility and processed-dataset guide

This public repository is the reproducibility release for the Radiation Shielding Research project. It contains the canonical Phase I and Phase II notebooks, the public result trees, the Phase I verification release, and a GitHub-safe byte-verifiable representation of the datasets needed to reconstruct the local scientific data tree.

This README is written for **external scientists who want to reproduce the published Phase I analysis, inspect the generated evidence, run the Phase II mathematical analysis, and use or understand the reconstructed datasets**. It is intentionally not a development log. Historical debugging, superseded source-construction experiments, private working material, and internal audit history are not part of this public release.

The public release contains both Phase I and Phase II material. The main methodological focus of this document remains Phase I because Phase I defines the frozen scientific corpus and physics baseline inherited by Phase II.

---

## 1. Project phases

### Phase I — external truth and conventional physics baseline

Phase I builds the scientific reference layer used by the rest of the project. It combines experimental, evaluated/reference, simulated, and derived datasets while preserving the evidence class of each source. It then performs conventional Geant4 calculations and source-model validation before constructing ordered response fields for downstream mathematical analysis.

The major Phase I tasks are:

1. assemble and normalize the external reference corpus;
2. construct a common evidence-aware schema without collapsing physically different observables into one meaning;
3. register the production photon source models;
4. validate the required 6, 10, 15, 16, and 18 MV photon sources against independent measurements;
5. reproduce the NIST ordinary-concrete photon benchmark;
6. reproduce the JAERI/TIARA 43 and 68 MeV neutron transmission and dose-equivalent benchmarks;
7. generate controlled shielding-response sweeps and Monte Carlo residual fields;
8. integrate the expanded photon, neutron, electron, photonuclear, and shielding-spectra corpus;
9. construct ordered reference and residual fields for Phase II;
10. apply the Phase I exit gate.

The current Phase I result is stored under `results/phase1/`, and the release artifacts are stored under `releases/phase1/`.

### Phase II — mathematical structure discovery

Phase II consumes the frozen Phase I corpus and ordered response fields rather than redefining the Phase I experimental evidence. The current public Phase II notebook is:

```text
notebooks/phase2/01_Phase2_RELEASE_v13.ipynb
```

Supporting public Phase II outputs are under:

- `results/phase2/`

The Phase II workflow evaluates mathematical structures against radiation-response fields. It is downstream of the Phase I physics baseline and does not replace the conventional validation described here. Private development notes and intermediate research workspaces are intentionally not part of the public release.

For `P005_SIM` (PSSD), Phase II reads the authoritative canonical table directly:

```text
data/processed/PSSD/pssd_canonical_spectra.parquet
```

The Phase II object is `P005_SIM.TRAJ_SPEC`: at fixed `(element, incident_energy_MeV)`, the ordered coordinate is `depth_MFP`, and each state is the complete outgoing-energy spectrum `relative_flux(outgoing_energy_MeV)`. The canonical PSSD therefore supplies a **16-depth × 30-energy** vector-valued trajectory for each element/incident-energy condition.

Phase II does not reconstruct PSSD axes from filenames, paths, wide `Flux_*` columns, or compatibility-table heuristics. It also does not invent additional shielding depths by interpolation or reinterpret outgoing-energy bins as the trajectory coordinate merely to satisfy estimator sample-count requirements. Methods whose declared minimum support exceeds the 16 available depth states are expected to return a conditional/unsupported result rather than a synthetic result.

### Later phases

Any later project phases are downstream of the frozen Phase I evidence and the Phase II mathematical-analysis work. They are outside the scope of this README unless a dedicated public reproduction procedure is added to the repository.

---

## 2. Public repository layout

The public release intentionally contains only the reproducibility surfaces needed by external users:

| Path | Purpose |
|---|---|
| `README.md` | Primary public scientific and reproducibility guide. |
| `LICENSE` | Repository license. |
| `environment.yml` | Portable Python/analysis environment specification. |
| `environment-lock-linux-64.txt` | Linux-64 lock record for the Python/analysis environment. |
| `environment-geant4.yml` | Native Geant4 environment specification. |
| `environment-geant4-lock-linux-64.txt` | Linux-64 lock record for the Geant4 environment. |
| `requirements-full.txt` | Full Python package inventory for reproducibility. |
| `datasets/github/` | GitHub-safe, byte-verifiable distribution representation of the scientific dataset. |
| `notebooks/phase1/01_Phase1.ipynb` | Canonical Phase I notebook. |
| `notebooks/phase2/01_Phase2_RELEASE_v13.ipynb` | Canonical Phase II notebook. |
| `results/phase1/` | Public Phase I results, Geant4 outputs, validation metrics, response fields, and final gate state. |
| `results/phase2/` | Public Phase II results. |
| `releases/phase1/` | Phase I verification/release archive. |

The authoritative scientific dataset is **not checked into Git as `data/`**. Instead, the repository stores the distributable representation under `datasets/github/`. The reconstruction helper and manifest in that tree are the authoritative mechanism for restoring and verifying the dataset locally.

Private development trees, historical archives, source-reconstruction workspaces, helper-script collections, and internal documentation directories are intentionally omitted from the public release. Their absence does not affect the canonical Phase I or Phase II notebooks or the published result trees.

---

## 3. Dataset installation: one authoritative instruction location

The public Git checkout does not contain a tracked `data/` directory. The scientific dataset is distributed in GitHub-safe form under:

```text
datasets/github/
```

Large files and very wide source directories may therefore appear in the GitHub tree as split objects or deterministic directory bundles rather than in their final local form.

The authoritative distribution/reconstruction instructions are maintained in the project's second maintained README:

```text
datasets/github/_distribution_meta/README.md
```

The reconstruction helper is:

```text
datasets/github/_distribution_meta/recombine_dataset.py
```

To reconstruct the GitHub distribution safely into the helper's default local-only `datasets/reconstructed/` tree:

```bash
python datasets/github/_distribution_meta/recombine_dataset.py
```

The exact reconstruction, verification, and local working-tree instructions are maintained in `datasets/github/_distribution_meta/README.md`. Follow that document before running analyses that require the reconstructed `data/...` paths shown throughout this guide.

If an authoritative local `data/` tree is already present, it can be verified without modifying it:

```bash
python datasets/github/_distribution_meta/recombine_dataset.py \
  --output data \
  --verify-only
```

`datasets/github/` is a distribution representation, not the analysis data root. Scientific code should use the reconstructed local dataset according to the distribution README.

The root `README.md` and `datasets/github/_distribution_meta/README.md` are the two project-maintained README instruction surfaces. No separate `DATASET_RECONSTRUCTION.md` is required.

---

## 4. Software and computational requirements

Phase I combines Python/Jupyter analysis with native multithreaded Geant4 transport.

### Python layer

The Phase I notebook uses standard-library Python together with:

- NumPy
- pandas
- PyArrow
- Matplotlib
- IPython/Jupyter

The canonical Phase I notebook is:

```text
notebooks/phase1/01_Phase1.ipynb
```

The canonical Phase II notebook is:

```text
notebooks/phase2/01_Phase2_RELEASE_v13.ipynb
```

Portable environment records are provided at the repository root:

```text
environment.yml
environment-lock-linux-64.txt
environment-geant4.yml
environment-geant4-lock-linux-64.txt
requirements-full.txt
```

The notebooks resolve the repository automatically when run from inside the project tree. A nonstandard checkout location can be supplied with the `RADIATION_SHIELDING_REPO` environment variable.

### Geant4 layer

The native Phase I implementation is C++17 and requires a **multithreaded Geant4 build**. The recorded Route-C campaign used Geant4 11.4.2 with multithreading enabled.

The experiment contract requires:

- C++17;
- Geant4 multithreading;
- exactly **16 Geant4 worker threads** for the frozen Phase I transport contracts;
- no serial fallback for required stochastic runs;
- at least **1,000,000 histories** for required stochastic Phase I runs;
- higher run-specific neutron statistics where specified;
- **100,000,000 histories per beam** for the paired 6/10/15/16/18 MV Route-C photon campaign.

The main native Geant4 source is:

```text
results/phase1/geant4_cpp/main.cpp
```

with build definition:

```text
results/phase1/geant4_cpp/CMakeLists.txt
```

A standard build is:

```bash
cmake -S results/phase1/geant4_cpp -B results/phase1/geant4_cpp/build
cmake --build results/phase1/geant4_cpp/build -j16
```

The resulting executable supports these scientific modes:

```text
phase1_geant4_runner info
phase1_geant4_runner p001 <target.csv> <out.csv>
phase1_geant4_runner neutron <config.ini>
phase1_geant4_runner photon_pdd <config.ini>
```

Prepared run configurations are stored under:

```text
results/phase1/geant4_templates/run_configs/
```

and photon-PDD configurations under:

```text
results/phase1/geant4_templates/run_configs/production_photon_pdd/
```

---

## 5. Reproducing the final Phase I analysis

Once the local `data/` tree has been reconstructed according to the authoritative data instructions, the final notebook can be executed in audit/evaluation mode using the existing transport outputs:

```bash
cd /path/to/Radiation-Shielding-Research

RADIATION_SHIELDING_REPO="$PWD" \
PHASE1_MODE=audit \
STEP2_ROUTE_C_MODE=audit \
python -m jupyter nbconvert \
  --to notebook \
  --execute \
  notebooks/phase1/01_Phase1.ipynb \
  --output phase1_audit_run.ipynb \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.allow_errors=False
```

This evaluates the corpus, validation results, conventional benchmarks, response fields, and Phase I exit gate without launching a new high-statistics transport campaign.

The principal completion files are:

```text
results/phase1/phase1_complete_status.json
results/phase1/phase1_exit_gate.csv
results/phase1/phase1_audit_summary.json
results/phase1/phase1_implementation_status.csv
```

A complete Phase I result requires every active row in:

```text
results/phase1/phase1_exit_gate.csv
```

to pass.

The verification archive is written under:

```text
releases/phase1/
```

---

## 6. Reproducing the Phase II analysis

After reconstructing the scientific dataset according to the distribution README, the canonical Phase II notebook can be executed from the repository root:

```bash
python -m jupyter nbconvert \
  --to notebook \
  --execute \
  notebooks/phase2/01_Phase2_RELEASE_v13.ipynb \
  --output phase2_run.ipynb \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.allow_errors=False
```

The public Phase II result tree is under:

```text
results/phase2/
```

The checked-in canonical notebook is intentionally stored without execution outputs or execution counts. The public result tree is the release record of the completed analysis.

---

## 7. Reproducing the Geant4 experiments

The final notebook is the analysis/orchestration record. The frozen C++ sources, run configurations, and output contracts under `results/phase1/` are the direct reproducibility artifacts for the transport experiments.

### 7.1 P001 — NIST ordinary-concrete photon attenuation

P001 compares deterministic Geant4 electromagnetic attenuation coefficients with the NIST ordinary-concrete reference.

Reference data:

```text
data/benchmarks/photon/nist_ordinary_concrete_attenuation.csv
data/benchmarks/photon/nist_ordinary_concrete_composition.csv
data/benchmarks/photon/nist_ordinary_concrete.npz
data/benchmarks/photon/metadata.json
```

The attenuation table preserves all 53 NIST rows, including duplicate energies associated with absorption-edge discontinuities. Derived linear attenuation coefficients use a concrete density of 2.300 g/cm³ and are convenience values rather than separately published NIST values.

Primary generated targets and comparisons are under:

```text
results/phase1/validation_targets/
results/phase1/geant4_raw/
results/phase1/conventional_baseline/
```

The frozen P001 acceptance criteria are:

- maximum absolute relative error ≤ 5%;
- median absolute relative error ≤ 2%.

Absorption-edge rows are handled using the repository's deterministic edge-probe procedure rather than by silently collapsing duplicate energies.

### 7.2 N001/N002 — JAERI/TIARA neutron transmission

The neutron benchmark reproduces the JAERI/TIARA concrete-shielding measurements for 43 and 68 MeV proton-produced neutron fields.

Reference data:

```text
data/benchmarks/neutron/jaeri_tiara_concrete_composition.csv
data/benchmarks/neutron/jaeri_tiara_experiment_setup.csv
data/benchmarks/neutron/jaeri_tiara_source_spectra.csv
data/benchmarks/neutron/jaeri_tiara_bc501a_transmission_spectra.csv
data/benchmarks/neutron/jaeri_tiara_dose_equivalent.csv
data/benchmarks/neutron/jaeri_tiara_concrete.npz
data/benchmarks/neutron/metadata.json
```

The BC501A transmission table contains the published spectrum points in long form by source energy, concrete thickness, detector position, and neutron-energy bin. Missing published off-axis bins are not invented or interpolated.

Important reproduction rules:

- the photon and neutron benchmark concrete compositions are different and must not be substituted for one another;
- measured neutron source spectra are preserved;
- absolute source normalization comes from the experiment-setup data;
- no free post-hoc normalization of neutron Monte Carlo spectra is permitted;
- spectral agreement is evaluated only after the Monte Carlo statistical-adequacy gate passes.

The neutron spectral-statistics adequacy requirements are:

- positive-bin coverage ≥ 90%;
- median relative Monte Carlo uncertainty ≤ 40%;
- 90th-percentile relative Monte Carlo uncertainty ≤ 75%.

Once that statistics gate is satisfied, the conventional neutron-spectrum agreement criterion is:

- median multiplicative disagreement factor ≤ 1.5;
- at least 80% of reference-positive bins within a factor of two.

The ICRP-21 dose-equivalent comparison uses the same factor limits:

- median multiplicative disagreement factor ≤ 1.5;
- at least 80% of compared values within a factor of two.

Prepared neutron and normalization configs are under:

```text
results/phase1/geant4_templates/run_configs/
```

Required neutron campaigns use at least the general 1,000,000-history floor, with run-specific high-statistics floors retained in the configuration and planning files. The Phase I policy uses 10,000,000 histories for moderate high-statistics cases and 100,000,000 histories for severe cases.

### 7.3 Controlled shielding-response sweep

A separate fixed-geometry concrete-thickness sweep is retained for mathematical-response analysis. It should not be confused with the historical JAERI benchmark geometries, which vary with thickness and collimation. The controlled sweep provides a consistent response field for Phase II mathematical analysis.

Its configuration and outputs are kept with the main Geant4 run configuration and result trees.

---

## 8. Production photon source models and Step 2C

### 8.1 Model scope

The required Phase I production sources are the frozen:

- 6 MV
- 10 MV
- 15 MV
- 16 MV
- 18 MV

photon-energy distributions.

The 3 MV source may exist in the simulator but is not a blocking Phase I production source.

These production sources are **one-dimensional photon-energy distributions**, not unique clinical x/y/energy/angle phase-space models. Validation is therefore aligned with the physical scope of the model being validated.

### 8.2 Current Step 2C operational validation rule

A required source passes operational energy-source validation when at least one independent measured, non-TVL modality appropriate to the 1-D source scope passes its predefined numerical criterion.

Qualifying modalities are:

1. **Measured PDD**, which tests beam-quality and depth-dose behavior; or
2. **Independent measured photon spectrum**, when both support and shape criteria pass.

A contradictory strict machine-matched spectral reference remains a veto.

Lateral dose profiles are retained as diagnostics of the factorized spatial surrogate. They are **not** used as a pass/fail test of the 1-D energy distribution because profile shape depends strongly on spatial and angular source structure that is not uniquely specified by a 1-D spectrum.

Cross-machine or different-field spectral evidence is retained as a comparability qualifier rather than silently treated as machine-identical evidence.

No production spectrum, acceptance threshold, smoothing rule, normalization rule, or measured dataset is modified in order to obtain a validation result.

### 8.3 Measured-spectrum criterion

Production source probabilities are photon-number sampling masses. Spectral validation converts them to energy-fluence mass using:

```text
energy_fluence_mass_i = energy_center_MeV_i × probability_mass_bin_i
```

The two-sided support requirement is:

- ≥ 95% of production energy fluence covered by the measured reference; and
- ≥ 95% of measured-reference energy fluence covered by positive production support.

The frozen spectral-shape thresholds are:

- Jensen-Shannon divergence ≤ 0.08 bits;
- total variation ≤ 0.25;
- cosine similarity ≥ 0.95.

### 8.4 Measured-PDD criterion

Each reference PDD and Monte Carlo PDD is normalized to its own dmax = 100%.

Evaluation is post-dmax and requires at least 10 post-dmax points.

The frozen PDD thresholds are:

- median absolute difference ≤ 2.5 percentage points;
- 90th-percentile absolute difference ≤ 5.0 percentage points;
- absolute PDD(10 cm) difference ≤ 3.0 percentage points.

### 8.5 Current qualifying evidence by source

The current Phase I Step 2C evaluation qualifies the required source family as follows:

| Source | Qualifying independent evidence |
|---|---|
| 6 MV | measured PDD |
| 10 MV | measured PDD |
| 15 MV | measured PDD |
| 16 MV | measured spectrum support + shape |
| 18 MV | measured PDD |

The current operational results are written to the source-validation area under:

```text
results/phase1/source_models/validation/
```

with the active Step 2C operational table and policy files generated by the notebook.

---

## 9. Route-C paired PDD/profile experiment

The paired Route-C campaign provides a second, multi-observable look at the factorized source surrogate.

The frozen campaign uses:

- 6, 10, 15, 16, and 18 MV production sources;
- a 40 × 40 cm² field;
- the same simulated events for PDD and the 10-cm profile;
- exactly 16 Geant4 worker threads;
- at least 100,000,000 histories per beam;
- a factorized water-surface source-plane model.

The frozen C++ runner is:

```text
results/phase1/source_models/validation/route_c_paired_transport_v1/step2_route_c_paired_runner.cpp
```

Its command-line interface is:

```text
step2_route_c_runner info
step2_route_c_runner paired <config.ini>
```

The five frozen beam configurations are under:

```text
results/phase1/source_models/validation/route_c_paired_transport_v1/configs/
```

After building the runner with a multithreaded Geant4 C++17 toolchain, each configuration can be executed with the `paired` mode. The exact compiler/toolchain record from the reference campaign is retained in the Route-C log directory.

### Profile diagnostic convention

Profiles are normalized at the published central axis, x = 0, to 100%.

The profile diagnostic uses:

- reference inclusion threshold ≥ 10% of CAX;
- plateau region `|x| ≤ 16 cm` with reference ≥ 80%;
- plateau 90th-percentile absolute difference threshold of 2.0 percentage points;
- bilateral 20%, 50%, and 80% edge-crossing distance-to-agreement threshold of 2.0 mm.

No mirroring, profile symmetrization, smoothing, or fitted lateral translation is permitted.

These profile metrics are reported to characterize the factorized spatial surrogate; they are not the Step 2C gate for the 1-D energy spectrum.

---

## 10. Step 2 multi-observable holdout datasets

The canonical processed holdouts are under:

```text
data/processed/production_source_validation/step2_holdouts/
```

The ten primary measured observables are:

| Energy | PDD reference | 10-cm profile reference |
|---|---|---|
| 6 MV | Ben Hdech, Figure 62 | Ben Hdech, Figure 62 bilateral profile |
| 10 MV | Araki 2005, Figure 3(d) | Araki 2005, Figure 5(d) |
| 15 MV | Tian et al. 2017, Figure 4 | Tian Figure 5(b2), 10-cm cross-line |
| 16 MV | Ben Hdech, Figure 62 | Ben Hdech, Figure 62 bilateral profile |
| 18 MV | Ben Hdech, Figure 62 | Ben Hdech, Figure 62 bilateral profile |

Processed files are:

```text
6MV_PDD.csv
6MV_PROFILE_10cm.csv
10MV_PDD.csv
10MV_PROFILE_10cm.csv
15MV_PDD.csv
15MV_PROFILE_10cm_CROSSLINE.csv
16MV_PDD.csv
16MV_PROFILE_10cm.csv
18MV_PDD.csv
18MV_PROFILE_10cm.csv
step2_holdouts_long.csv
step2_holdout_dataset_registry.csv
step2_holdout_intrinsic_qc.csv
```

The raw source evidence used to create them is under:

```text
data/raw/production_source_validation/step2_holdouts/
```

including source documents, exact graphical objects, source freezes, and provenance material.

### Digitization rules

The holdouts were frozen before project predictions were viewed. The extraction procedure does not permit:

- smoothing;
- symmetrization;
- profile mirroring;
- fitted horizontal translation;
- post-hoc acceptance-threshold changes.

PDDs are normalized to their own dmax = 100%.

Profiles are normalized to the published central-axis value at x = 0 = 100%.

### Source-specific extraction details

**Ben Hdech — 6, 16, and 18 MV**

Figure 62 is native vector artwork. The black curves are experimental data; red curves are the publication's Monte Carlo data. The experimental vector paths used for the 40 × 40 cm² holdouts are:

- 6 MV PDD: path 255
- 6 MV profile: path 467
- 16 MV PDD: paths 716 + 717 + 718 + 719
- 16 MV profile: paths 942 + 943
- 18 MV PDD: path 1187
- 18 MV profile: path 1399

The four 16 MV PDD path fragments join exactly in rendered vector coordinates. The two 16 MV profile fragments are preserved as published and are not snapped together.

Ben Hdech Figure 62 uses display factors by field size. For the 40 × 40 cm² curves the display factor is `k = 0.4`:

- for the PDD, the factor is removed before final dmax normalization;
- for the profile, the multiplicative factor cancels under CAX normalization.

**Araki — 10 MV**

The exact PDF-embedded raster objects are used:

- Figure 3(d): 10 MV, 40 × 40 cm² PDD;
- Figure 5(d): 10 MV, 40 × 40 cm² profile at 10 cm depth.

Only the measured solid curve is treated as holdout data.

**Tian — 15 MV**

The exact PDF-embedded JPEG objects are used:

- Figure 4: 15 MV FF PDD;
- Figure 5(b2): 10-cm cross-line profile.

The measurement is the red solid curve; the 40 × 40 cm² curve is the outermost/widest profile in the selected panel.

The shallow Tian profile is intentionally excluded because the publication gives inconsistent depth metadata: Figure 5 / Section 3.A.2 reports 2.8 cm, while Table 4 reports 2.1 cm for the 40 × 40 cm² field. No reconciliation is assumed.

---

## 11. Production-source validation corpus

The broader production-source validation corpus is under:

```text
data/processed/production_source_validation/
```

It contains:

- measured 40 × 40 cm² photon PDD validation data for 6, 10, 15, 16, and 18 MV;
- a measured Mathew 15 MV neutron-spectrum dataset;
- a measured Zamorano thermal-neutron-fluence dataset;
- a Zamorano/Digital.CSIC PHITS head-neutron spectrum retained explicitly as a **simulated code-to-code reference**, not as an experiment.

The principal validation dataset identifiers represented in the processed registry include:

```text
PVAL_6MV_TRUEBEAM_40x40_LANDAUER
PVAL_10MV_TRUEBEAM_40x40_LANDAUER
PVAL_15MV_TRUEBEAM_40x40_LANDAUER
PVAL_16MV_CLINAC21IX_40x40_PACYNIAK
PVAL_18MV_ONCOR_40x40_SAWKEY
NVAL_15MV_TRUEBEAM_STX_MATHEW
NVAL_15MV_ZAMORANO_THERMAL
NVAL_15MV_ZAMORANO_PHITS_SOURCE
```

These names are scientific dataset identifiers, not top-level directories in the public Git checkout. Use the reconstructed processed-data registry to resolve their canonical files.

The combined photon-PDD table is:

```text
all_photon_pdd_validation_6_10_15_16_18MV.csv
```

The dataset registry and processing quality-control files in this directory are the preferred entry points for programmatic use.

The 16 MV Pacyniak PDD retains its source-reporting geometry caveat. Where the exact experimental SSD is unresolved, the comparison should be interpreted as normalized shape evidence rather than silently upgraded to an exact-geometry claim.

---

## 12. Expanded processed corpus

Phase I also registers an expanded corpus under `data/processed/`. The evidence class of each dataset is scientifically important and must be preserved.

| Dataset family | Phase I role / evidence class |
|---|---|
| `BROAD_BEAM_PHOTON` | experimental broad-beam photon shielding/attenuation data |
| `IAEA_EXFOR` | experimental photonuclear reaction measurements |
| `IAEA_PD2019` | evaluated photonuclear reaction data |
| `NIST_ESTAR` | evaluated electron transport/stopping-power data |
| `PSSD` | independent simulated photon shielding spectra in an authoritative canonical long-form representation; not experimental truth |
| `RB2000164` | experimental ISIS neutron-transmission data |
| `RB2000209` | experimental-derived standard-monitor transmission proxy |
| `production_source_validation` | measured and simulated source-validation data described above |

### Evidence semantics

The project does not silently treat experimental, evaluated, simulated, and derived/proxy evidence as equivalent.

In particular:

- PSSD is a simulated independent reference and is not experimental truth;
- PD-2019 is a photonuclear reaction-data layer, not an XCOM-like photon attenuation table;
- EXFOR observations preserve target, reaction channel, experiment identity, uncertainty, and dependence/status metadata;
- correlated or dependent EXFOR observations are not counted as independent measurements without a defensible statistical treatment;
- missing experimental uncertainty is kept missing rather than invented to create a z-score;
- RB2000209 is treated as transmission data only; a removal cross section is not claimed unless the required thickness definition is physically grounded;
- the current RB2000209 archive product is not labeled as a GEM product;
- canonical strict long-form experimental data are preferred over convenience matrices created by resampling.

### PSSD canonical representation

The authoritative processed PSSD spectral product is:

```text
data/processed/PSSD/pssd_canonical_spectra.parquet
```

It is a canonical long-form numerical table with exactly these fields:

```text
particle
atomic_number
element
incident_energy_MeV
depth_MFP
outgoing_energy_MeV
relative_flux
```

The current canonical table contains:

- **971,520 rows**;
- **92 elements**;
- **22 incident photon energies**;
- **16 shielding-depth states** in `depth_MFP`;
- **30 outgoing-energy bins**;
- **2,024 element × incident-energy trajectories**.

For each fixed `(element, incident_energy_MeV)` condition, the 480 rows form one complete **16 × 30 depth-by-outgoing-energy spectral tensor**. `atomic_number` is retained as canonical element identity metadata.

The preserved source-derived wide table is:

```text
data/processed/PSSD/tables/PSSD_Usable_Data_relative_factors_data.csv.parquet
```

Its `Flux_*MeV` columns were normalized upstream when the canonical table was built; Phase II does not repeat that decoding. The original source-derived table remains preserved rather than being rewritten into an analysis-specific schema.

For compatibility with older code, the project also retains:

```text
data/processed/PSSD/pssd_semantic_spectra.parquet
```

That file is a **compatibility product**, not the authoritative PSSD endpoint. New downstream analysis should use `pssd_canonical_spectra.parquet`.

### Ordered-field distinction

For downstream analysis:

- RB2000164 has an exact common experimental energy grid;
- RB2000209 strict long-form observations do not share an identical per-sample energy grid, so any matrix representation is explicitly a derived/resampled convenience product;
- PSSD is represented canonically as depth-ordered complete spectra: each `(element, incident_energy_MeV)` condition defines a 16-state `depth_MFP` trajectory whose state vector contains 30 outgoing-energy `relative_flux` values;
- PSSD pure-element spectra are not composition-weighted to manufacture synthetic concrete spectra, and missing depth states are not created merely to satisfy downstream estimator support requirements;
- PD-2019 is indexed as isotope/reaction/energy curves.

---

## 13. Common scientific schema

The Phase I common schema is designed to make heterogeneous radiation data queryable while preserving the physical meaning of each observable.

Conceptually, records retain:

```text
particle
incident energy
outgoing/scored energy
spatial or shielding coordinate
material
observable
source definition
uncertainty
scientific evidence class
```

The schema is not intended to imply that attenuation coefficients, neutron transmission, stopping powers, dose, PDD, photonuclear cross sections, and simulated spectra are dimensionally interchangeable. The observable and evidence fields must remain part of any downstream analysis.

Common-schema outputs are under:

```text
results/phase1/common_schema/
```

Ordered physical response fields are under:

```text
results/phase1/ordered_response_fields/
```

---

## 14. Using the benchmark datasets directly

### Photon benchmark example

```python
import pandas as pd

photon = pd.read_csv(
    "data/benchmarks/photon/nist_ordinary_concrete_attenuation.csv"
)
print(photon.head())
```

### Neutron benchmark example

```python
import pandas as pd

transmission = pd.read_csv(
    "data/benchmarks/neutron/jaeri_tiara_bc501a_transmission_spectra.csv"
)
source = pd.read_csv(
    "data/benchmarks/neutron/jaeri_tiara_source_spectra.csv"
)
setup = pd.read_csv(
    "data/benchmarks/neutron/jaeri_tiara_experiment_setup.csv"
)
```

### Combined NumPy benchmark archive

```python
import numpy as np

data = np.load("data/benchmarks/radiation_benchmarks_core.npz")
print(data.files)
```

Important: normalized neutron source spectral shapes and absolute source flux are separate pieces of experimental information. Use the experiment-setup table for the corresponding absolute source normalization.

---

## 15. Using the processed datasets programmatically

The processed-data registry files should be preferred over hard-coding assumptions about file names or dataset semantics.

For the expanded corpus, the notebook reads the processed dataset registry and then loads canonical Parquet files by dataset identifier. A minimal pattern is:

```python
import pandas as pd

registry = pd.read_parquet(
    "data/processed/shielding_dataset_registry.parquet"
)
print(registry[["dataset_id", "canonical_file", "evidence_type"]])
```

For the production-source validation corpus, use:

```text
data/processed/production_source_validation/shielding_validation_dataset_registry.csv
```

and for Step 2 holdouts use:

```text
data/processed/production_source_validation/step2_holdouts/step2_holdout_dataset_registry.csv
```

The dataset-specific metadata files should be read together with the numerical tables when interpreting machine, beam, geometry, normalization, or evidence class.

---

## 16. Scientific rules that must be preserved in independent reproduction

To reproduce the experiments as performed here, do not change the following rules after looking at the comparison results:

- do not mutate the production photon spectra to improve validation agreement;
- do not tune numerical acceptance thresholds after viewing the holdout results;
- do not smooth measured holdout curves;
- do not symmetrize or mirror measured profiles;
- do not fit a lateral translation to improve a profile comparison;
- do not invent missing experimental bins or uncertainties;
- do not apply free post-hoc normalization to neutron spectra;
- do not treat simulated reference datasets as experimental measurements;
- do not merge distinct reaction channels or nuclides in EXFOR/PD-2019 comparisons;
- do not treat a 1-D photon energy spectrum as if it uniquely specified a complete clinical phase space;
- preserve the source-specific material composition and geometry for benchmark transport;
- preserve the required history counts and 16-worker multithreaded transport contracts.

These constraints are part of the scientific method of the Phase I reproduction, not formatting preferences.

---

## 17. Where to inspect the final results

For an external reviewer, the most useful Phase I result files are:

```text
results/phase1/phase1_complete_status.json
results/phase1/phase1_exit_gate.csv
results/phase1/phase1_audit_summary.json
results/phase1/phase1_implementation_status.csv
results/phase1/conventional_baseline/phase1_conventional_baseline_summary.csv
results/phase1/conventional_baseline/neutron_spectral_statistics_adequacy.csv
results/phase1/source_models/validation/phase1_source_validation_case_results.csv
results/phase1/source_models/validation/step2c_operational_validation_v2.csv
results/phase1/source_models/validation/route_c_paired_transport_v1/step2_route_c_pdd_metrics.csv
results/phase1/source_models/validation/route_c_paired_transport_v1/step2_route_c_profile_metrics.csv
results/phase1/ordered_response_fields/
```

The corresponding raw Geant4 outputs and provenance records are under:

```text
results/phase1/geant4_raw/
```

The final portable Phase I verification archive is under:

```text
releases/phase1/
```

---

## 18. Raw evidence and reproducibility provenance

Raw source documents, machine-readable reference material, and the public provenance needed to trace the processed datasets are distributed through `datasets/github/` and become available through the reconstructed scientific data tree.

For the Step 2 photon holdouts, the reconstructed raw evidence is under:

```text
data/raw/production_source_validation/step2_holdouts/
```

The corresponding processed holdouts and registries are under:

```text
data/processed/production_source_validation/step2_holdouts/
```

Public Phase I provenance and validation outputs are retained under `results/phase1/`, including the source-validation records, Geant4 outputs, ordered response fields, manifests, and verification summaries needed to inspect the released analysis.

Private development notes, historical debugging packages, and source-reconstruction workspaces are intentionally not published as part of this release. The public README therefore documents the scientifically relevant extraction rules and evidence semantics directly rather than referring users to private working directories.

Source-specific license, citation, metadata, or provenance files contained in the reconstructed dataset should be consulted when interpreting or redistributing third-party source material.

---

## 19. Expected Phase I completion state

A successful reproduction of the current Phase I analysis should produce:

- a passing `phase1_exit_gate.csv` for all active Phase I requirements;
- `phase1_complete = true` in the completion status;
- operational Step 2C qualification for all five required photon sources;
- conventional P001, N001, and N002 benchmark outputs and agreement summaries;
- ordered reference and Monte Carlo residual fields;
- a Phase I verification archive under `releases/phase1/`.

The scientific interpretation should always distinguish what has been validated:

- the production photon sources are validated as **1-D energy-source models** using measured observables appropriate to that scope;
- the Route-C lateral profiles characterize the **factorized spatial surrogate**;
- neither of those claims should be inflated into a statement that a unique full clinical accelerator phase space has been reconstructed.

---

## 20. Consolidated documentation policy

This root README is the primary public-facing scientific reproduction guide for this release.

The project maintains two README instruction surfaces:

1. `README.md` — consolidated scientific context, experiment, dataset, validation, Phase II interface, and interpretation guidance;
2. `datasets/github/_distribution_meta/README.md` — GitHub dataset-distribution, reconstruction, verification, and upload-layout instructions.

No separate project-maintained dataset-reconstruction README is required. Third-party source archives may still contain upstream README or license files; those files document their original sources and are not additional project workflow instructions.

