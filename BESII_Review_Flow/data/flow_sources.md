# Flow Excitation-Function Data Sources

Compiled on 2026-04-29 for `dv1/dy` and `v2` vs. collision energy.

## Output Tables

- `flow_dv1dy_vs_energy.csv`: midrapidity directed-flow slope `dv1/dy` vs. `sqrt(s_NN)`.
- `flow_dv1dy_low_energy_expanded.csv`: low-energy `dv1/dy` additions from HADES/SIS18 and STAR fixed-target HEPData tables.
- `flow_v2_cross_experiment_vs_energy.csv`: historical cross-experiment near-midrapidity `v2` excitation function, vector-digitized from the NA49 publication EPS source.
- `flow_v2_star_bes_eta0_vs_energy.csv`: STAR charged-hadron `v2{EP}` near `eta=0`, 10-40% centrality, from STAR public data tables.
- `flow_v2_low_energy_expanded.csv`: expanded `v2` excitation function with SIS/AGS/FXT/SPS/RHIC comparison points from the STAR fixed-target HEPData compilation.

## Primary Sources

- STAR directed flow:
  - L. Adamczyk et al. (STAR), "Beam-Energy Dependence of Directed Flow of Protons, Antiprotons and Pions in Au+Au Collisions", Phys. Rev. Lett. 112, 162301 (2014), DOI: https://doi.org/10.1103/PhysRevLett.112.162301
  - HEPData record: https://doi.org/10.17182/hepdata.105867
  - Tables used: HEPData Table 5 (`p`, `pbar`, `pi+`, `pi-`) and Table 6 (`net-proton`).

- E895 directed flow comparison points:
  - H. Liu et al. (E895), "Sideward Flow in Au+Au Collisions Between 2A GeV and 8A GeV", Phys. Rev. Lett. 84, 5488 (2000), DOI: https://doi.org/10.1103/PhysRevLett.84.5488
  - Values in the CSV were vector-digitized from `fig3_new.eps` in STAR arXiv:1401.3043, which plots E895 comparison points with the STAR excitation function.

- NA49 directed and elliptic flow:
  - C. Alt et al. (NA49), "Directed and elliptic flow of charged pions and protons in Pb+Pb collisions at 40A and 158A GeV", Phys. Rev. C 68, 034903 (2003), DOI: https://doi.org/10.1103/PhysRevC.68.034903
  - Directed-flow comparison points were vector-digitized from `fig3_new.eps` in STAR arXiv:1401.3043.
  - Historical `v2` excitation-function points were vector-digitized from `v2edep.eps` in NA49 arXiv:nucl-ex/0303001.

- STAR inclusive charged-hadron elliptic flow:
  - L. Adamczyk et al. (STAR), "Inclusive charged hadron elliptic flow in Au+Au collisions at sqrt(s_NN)=7.7-39 GeV", Phys. Rev. C 86, 054908 (2012), DOI: https://doi.org/10.1103/PhysRevC.86.054908
  - STAR public data page: https://drupal.star.bnl.gov/STAR/files/starpublications/190/data.html
  - STAR 200 GeV reference:
    - B. I. Abelev et al. (STAR), Phys. Rev. C 77, 054901 (2008), DOI: https://doi.org/10.1103/PhysRevC.77.054901
    - STAR public data page: https://drupal.star.bnl.gov/STAR/files/starpublications/108/data.html

- Historical `v2` comparison sources embedded in NA49 Fig. 24:
  - E877: J. Barrette et al., Phys. Rev. C 56, 3254 (1997), DOI: https://doi.org/10.1103/PhysRevC.56.3254
  - CERES/NA45: J. Slivova, Nucl. Phys. A 715, 615c (2003), DOI: https://doi.org/10.1016/S0375-9474(02)01545-2
  - STAR 130 GeV: C. Adler et al., Phys. Rev. C 66, 034904 (2002), DOI: https://doi.org/10.1103/PhysRevC.66.034904
  - STAR 200 GeV QM02 summary: L. Ray, Nucl. Phys. A 715, 45c (2003), DOI: https://doi.org/10.1016/S0375-9474(02)01412-4
  - PHENIX: S. Esumi et al., Nucl. Phys. A 715, 599c (2003), DOI: https://doi.org/10.1016/S0375-9474(02)01541-5
  - PHOBOS: S. Manly et al., Nucl. Phys. A 715, 611c (2003), DOI: https://doi.org/10.1016/S0375-9474(02)01544-0

- Low-energy and fixed-target additions:
  - HADES/SIS18 proton flow:
    - J. Adamczewski-Musch et al. (HADES), Eur. Phys. J. A 59, 80 (2023), DOI: https://doi.org/10.1140/epja/s10050-023-00936-6
    - HEPData record: https://doi.org/10.17182/hepdata.152804
  - STAR 3 GeV identified-particle flow:
    - M. S. Abdallah et al. (STAR), Phys. Lett. B 827, 137003 (2022), DOI: https://doi.org/10.1016/j.physletb.2022.137003
    - HEPData record: https://doi.org/10.17182/hepdata.110656
    - Tables used: Fig. 5a proton `dv1/dy` and Fig. 5b proton/FOPI `v2` comparison tables.
  - STAR fixed-target 4.5 GeV and the expanded `v2` excitation-function compilation:
    - M. S. Abdallah et al. (STAR), Phys. Rev. C 103, 034908 (2021), DOI: https://doi.org/10.1103/PhysRevC.103.034908
    - HEPData record: https://doi.org/10.17182/hepdata.95903
    - Tables used: Figure 11 (`dv1/dy`) and Figure 14 (`v2` excitation function).
  - FOPI SIS elliptic-flow source represented in the STAR Fig. 14 HEPData compilation:
    - A. Andronic et al. (FOPI), Phys. Lett. B 612, 173 (2005), DOI: https://doi.org/10.1016/j.physletb.2005.02.060
  - E895 AGS elliptic-flow source represented in the STAR Fig. 14 HEPData compilation:
    - C. Pinkenburg et al. (E895), Phys. Rev. Lett. 83, 1295 (1999), DOI: https://doi.org/10.1103/PhysRevLett.83.1295

## Caveats

- STAR `dv1/dy` rows from HEPData are direct tabulated values with published statistical and systematic uncertainties.
- E895 and NA49 `dv1/dy` comparison rows are vector-digitized from the STAR EPS source. The STAR caption notes that the older measurements have comparable but not identical cuts.
- Low-energy `dv1/dy` additions are intentionally kept in a separate table because HADES/SIS18 and STAR fixed-target centralities/acceptances are not identical to the STAR BES-I 10-40% table.
- The historical cross-experiment `v2` table is intended as a literature-excitation-function reproduction of NA49 Fig. 24. It mixes pion and charged-hadron measurements and uses near-midrapidity, mid-central selections that are similar but not identical across experiments.
- The expanded `v2` table is direct HEPData from the STAR fixed-target Fig. 14 compilation. It is broader than the first-pass NA49-based table and includes FOPI SIS, EOS/E895/E877 AGS/Bevalac, STAR FXT, SPS, RHIC points. It should be treated as an excitation-function comparison, not a uniform-analysis data set.
- The STAR BES `v2` table is cleaner internally: charged hadrons, 10-40% centrality, near `eta=0`, mostly `0.2 < pT < 2.0 GeV/c`; the 200 GeV row from the older STAR public table uses `0.15 < pT < 2.0 GeV/c`.

## Rebuild

Run:

```bash
python3 scripts/build_flow_literature_data.py
```
