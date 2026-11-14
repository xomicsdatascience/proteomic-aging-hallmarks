# A03 Program Annotation & Evidence Decisions

1. **No 1:1 Hallmark Claims**: Programs are named by what they actually measure in HepG2 whole-cell lysates (e.g. `PROG_DNA_REPAIR`, `PROG_OXPHOS`, `PROG_INNATE_TLR`), not as direct readouts of the 12 organismal hallmarks.
2. **Unmeasurable Hallmarks in HepG2 Monoculture**:
   - **Dysbiosis** cannot be measured in cultured HepG2 cells; `PROG_INNATE_TLR` is explicitly labeled as host innate TLR sensing.
   - **Stem cell exhaustion** cannot be assessed in transformed HepG2 hepatoma cells and was removed from the primary aging-program panel.
   - **Secreted SASP / Intercellular communication** are only weakly proxied in whole-cell lysates (`PROG_INFLAMMATORY`, `PROG_CELL_CELL_COMM`).
3. **Transcript-Derived Priors**: `DEMAGALHAES_AGING_UP` and `DEMAGALHAES_AGING_DN` are retained as RNA-derived aging priors (`PRIOR_DEMAGALHAES_UP`, `PRIOR_DEMAGALHAES_DN`), not proteomic aging clocks.
4. **Nuisance/Shared-Response Modules**: `NUIS_E2F_PROLIF`, `NUIS_G2M_PROLIF`, and `NUIS_APOPTOSIS` are frozen in A03 for A06 specificity/growth sensitivity analyses.
