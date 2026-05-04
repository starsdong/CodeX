# 3 GeV phi/K- ratio study

This directory packages the 3 GeV Au+Au phi/K- analysis artifacts without
including the original UrQMD source tree or compiled build products.

Contents:

- `data/phi_k_npart_events_3gev_10k.tsv`: event-level UrQMD 4.0 summary for
  the 10k-event test sample.
- `data/phi_k_centrality3_besii_matched_vs_npart_3gev_10k.tsv`: UrQMD
  centrality-binned phi/K- ratios for 0-10%, 10-40%, and 40-60%.
- `data/besii_measured_phi_k_ratio_vs_npart_3gev.tsv`: STAR BESII measured
  midrapidity dN/dy phi/K- ratios derived from HEPData ins1897327 Figure 3(a,b).
- `data/thermus_phi_k_ratio_rc_scan_3gev.tsv`: THERMUS SCE phi/K- ratios for
  Rc = 2, 3, 4, and 6 fm using the BESII_Review 3 GeV fit parameters.
- `data/hepdata/`: the HEPData CSV tables used for the measured K- and phi
  points.
- `data/thermus/`: compact THERMUS manifests, summaries, and point predictions
  for the Rc scan.
- `figures/`: final PNG/PDF comparison plot.
- `scripts/run_phi_k_npart_10k.sh`: runner/parser used to build the UrQMD
  event summary.
- `scripts/make_phi_k_ratio_plot.py`: regenerates the final comparison plot
  from the packaged TSV tables.

THERMUS parameters used for the horizontal lines:

- Ensemble: SCanonical
- T = 81.662316 MeV
- muB = 751.420185 MeV
- gammaS = 1.0
- R = 8.076258 fm
- dV/dy = 2206.577381 fm^3
- B/(2Q) charge constraint = 1.246835
- Rc scan: 2, 3, 4, and 6 fm

