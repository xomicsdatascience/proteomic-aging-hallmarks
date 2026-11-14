# A01 Design & Provenance Audit

- **Input checksums**:
  - `data/raw (1).h5ad`: SHA256 `ff3a156ab342488a60a248613e7dfdf8fd53689e0b5996df5af5f1aa30c8ad05` (`1232` samples x `8630` features)
  - `data/drug_details (1).xlsx`: SHA256 `099e5d6b3df7cf73b1bf193363b290bcc33ae7c2a6bb38f0567ce78a811fc5cb` (`173` rows, `173` unique `ProductName`, `166` unique `BaseName`)
- **Sample Reconciliation**:
  - Total observations: `1232` across 13 plates (`Blue1–Blue6`, `Red1–Red6`, `RedMix`).
  - `library_compound`: `1,029` samples (`172` unique compounds: `160` in Blue/Red plates with `n=5–6` across 6 plates; `12` in `RedMix` with `n=6` on 1 plate).
  - `positive_assay_control` (`DFO`): `24` samples (`2` per plate on `Blue1–6` and `Red1–6`; `0` on `RedMix`).
  - `vehicle_control` (`DMSO`): `23` samples (`2` per plate on `Blue1–6` and `Red1–6`, except `Blue3` has `1`; `0` on `RedMix`).
  - `analytical_qc` (`QC_sample`): `155` samples (`12` per plate on `Blue1–6` and `Red1–6`; `11` on `RedMix`).
  - `unlabeled_excluded`: `1` sample (`Blue5_F12`, `drug == 'Drug'`).
- **Feature Reconciliation**:
  - `raw (1).h5ad` contains `8,630` unique gene-symbol-indexed protein features (vs. `8,703` in parent manuscript).
  - `7,871` features have `<= 50%` zeros across all samples (`7,296` have `<= 20%` zeros). Zeros represent missing/below-detection values (`10.60%` overall).
- **Contrast Estimability**:
  - **Directly estimable (`161` contrasts)**: `160` Blue/Red library compounds + `DFO`, all replicated across 5–6 plates with matched within-plate `DMSO` controls.
  - **Reference-anchored only (`12` contrasts)**: `12` `RedMix` compounds confined to a single plate with `0` within-plate `DMSO` controls. Excluded from primary inference and analyzed only in a named QC-anchored sensitivity analysis.
