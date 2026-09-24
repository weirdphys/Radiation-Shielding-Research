PHASE I VERIFICATION PACKAGE
============================
Notebook revision: v12_18_0_step2c_operational_validation_amendment
Execution mode: audit
Full Phase I complete: True
Phase-II frozen-corpus handoff ready: True
Generated UTC: 2026-09-24T15:05:22.665576+00:00

Scientific execution policy:
- Native implementation: C++17
- Geant4 worker threads: exactly 16
- Minimum histories per stochastic run: 1000000
- Run-specific neutron transmission high-stat floors (retained unchanged in v12.12): moderate=10000000, severe=100000000
- Neutron spectral-statistics gate: positive-bin coverage >= 90%, median relative MC sigma <= 40%, p90 <= 75%
- Serial fallback: prohibited
- Neutron free post-hoc normalization: prohibited
- TVL/HVL: secondary diagnostic, not sole production-source validation
- Clinical PDD from a factorized 1-D surface spectrum is diagnostic only unless spatial/energy/angular phase-space is sufficiently specified.
- Production spectra count as exact only when both probability masses and the production energy grid are recovered.
- Measured-spectrum validation uses energy-fluence on both sides: >= 95% of production energy fluence must be covered by the reference AND >= 95% of reference energy fluence must be covered by positive production support.
- Production probability_mass_bin is photon-number sampling mass; Step 2C converts it to E*p energy-fluence mass before shape comparison.
- Cross-machine/unresolved-machine spectra remain explicit comparability qualifiers and do not by themselves promote the historical strict SPECTRAL_VALIDATED classification.
- Step 2C-v2 operational validation is model-scope aligned: each required 1-D energy source must pass at least one independent measured non-TVL modality appropriate to that scope.
- Lateral-profile disagreement is retained as a spatial-surrogate diagnostic and is non-gating for validation of the 1-D energy source.
- The original strict spectral/provenance Step 2C result is preserved separately as a historical non-pass; v12.18 does not rewrite it retroactively.
- Frozen JS/TV/cosine thresholds are unchanged from v12.14.
- Cross-revision PDD reuse is accepted only after source/reference/result hashes, PDD C++ semantics, scoring/config semantics, histories and output structure are revalidated; notebook-revision equality is not required.
- v12.17 keeps production-source mutation forbidden and the legacy PDD audit fail-closed, while STEP2_ROUTE_C_MODE=run permits only the frozen paired 100M-history Route-C PDD/profile campaign.
- Simulated construction-sanity spectra are permanently nonqualifying for Step 2C; independent measured holdouts remain reserved validation evidence.

Start audit with:
- phase1_complete_status.json
- phase1_exit_gate.csv
- phase1_audit_summary.json
- phase1_audit_checks.csv
- phase1_implementation_status.csv
- phase1_remaining_requirements.csv
- phase1_verification_bundle_summary.json
- phase1_artifact_manifest.json
- phase1_plot_manifest.json
- source_models/source_library_snapshot/production_source_library_snapshot_manifest.json
- geant4_raw/phase1_run_session.json (RUN mode)
- geant4_raw/phase1_geant4_toolchain.json (toolchain discovery/build evidence)
- geant4_raw/phase1_geant4_dataset_preflight.json (dataset discovery/install/verification evidence)
- geant4_logs/geant4_install_datasets.log (when installation was attempted)
- geant4_raw/geant4_run_manifest.csv (completed real runs)
- geant4_raw/provenance/ (completed real runs)
- conventional_baseline/source_normalization_gate.csv
- conventional_baseline/phase1_conventional_baseline_summary.csv
- conventional_baseline/neutron_spectral_statistics_policy.json
- conventional_baseline/neutron_spectral_statistics_adequacy.csv
- validation_targets/P001_NIST_ORDINARY_CONCRETE_edge_probe_policy.json
- conventional_baseline/P001_v12_16_1_audit_rehydration.json
- validation_targets/P001_NIST_ORDINARY_CONCRETE_edge_scan_targets.csv
- geant4_raw/P001_NIST_ORDINARY_CONCRETE_edge_scan_geant4.csv
- conventional_baseline/P001_NIST_ORDINARY_CONCRETE_edge_scan_audit.csv
- geant4_templates/v12_12_rerun_policy.json
- geant4_templates/v12_10_neutron_high_stat_required_history_plan.csv
- geant4_templates/v12_10_jaeri_sinbad_tally_rerun_plan.csv
- geant4_templates/v12_10_jaeri_source_phase_space_contract.json
- geant4_templates/v12_10_jaeri_sinbad_scoring_contract.json
- source_models/validation/v12_validation_corpus_integrity.csv
- source_models/validation/v12_validation_dataset_registry.csv
- source_models/validation/v12_photon_pdd_run_manifest.csv
- source_models/validation/v12_16_1_photon_pdd_reuse_audit.csv
- source_models/construction_audit/source_evidence_layer_contract.csv
- source_models/construction_audit/source_construction_contract.csv
- source_models/construction_audit/source_dependency_graph.csv
- source_models/construction_audit/source_structural_qa.csv
- source_models/construction_audit/public_sanity_reference_registry.csv
- source_models/construction_audit/public_mc_sanity_comparison.csv
- source_models/construction_audit/derived_source_reproducibility.csv
- source_models/construction_audit/source_reconstruction_decision.csv
- source_models/construction_audit/current_16MV_measured_shape_baseline.json
- source_models/construction_audit/v12_16_source_mutation_guard.json
- source_models/construction_audit/v12_16_source_construction_audit_summary.json
- source_models/validation/v12_5_photon_pdd_theory_contract.json
- source_models/validation/phase1_source_validation_case_results.csv
- source_models/validation/route_c_paired_transport_v1/step2_route_c_validation_summary.json
- source_models/validation/route_c_paired_transport_v1/step2_route_c_beam_results.csv
- source_models/validation/route_c_paired_transport_v1/step2_route_c_paired_run_manifest.csv
- source_models/validation/step2c_spectral_semantics_and_comparability_registry.csv
- source_models/validation/step2c_spectral_diagnostics_summary.csv
- source_models/validation/spectral_diagnostics/ (per-case energy-fluence comparison tables)
- geant4_templates/v12_15_step2c_semantics_support_comparability_policy.json
- documentation/<Phase-I scientific guide ODT> (when a companion guide is present at bundle creation)

A FULL Phase-I PASS is accepted only when every row in phase1_exit_gate.csv passes.
The original strict Step-2C spectral/provenance result remains available as historical evidence; the active Phase-I source gate uses the transparent Step-2C-v2 operational amendment recorded in source_models/validation/step2c_operational_validation_v2_policy.json.

The archive's own SHA-256 is written next to the ZIP after creation as
Phase1_v12_18_0_verification_bundle_sha256.json. That external hash file is
not embedded in the ZIP because doing so would create a self-referential archive hash.
