STEP 2C PUBLIC NUMERIC SPECTRUM INPUTS
=====================================

A publication/figure citation alone never passes Step 2C.
To activate an automatic measured_source_spectrum case:
1. Place a numerical CSV in public_numeric_references/.
2. Preserve the published/digitized energy coordinate in MeV and the spectral density.
3. Record the exact CSV SHA-256 in step2c_public_evidence_registry.csv.
4. Set machine_readable_numeric=True and numeric_reference_file to the CSV filename.
5. Keep used_to_construct_project_model=False and independent_of_project_model_fit=True.
6. Run the notebook; only a matching file/hash is appended to the Step 2C evaluator. v12.16 retains explicit quantity conversion, two-sided support checks, and comparability classification before any strict spectral PASS is possible.

Required default columns:
  energy_MeV
  normalized_energy_fluence_per_MeV

If a publication uses different columns or probability-mass semantics, declare those explicitly in the registry rather than silently converting them.