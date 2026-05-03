# Lambda Global Polarization Data

This folder collects published global-polarization measurements for Lambda and anti-Lambda hyperons in heavy-ion collisions. The main plot uses only `particle=Lambda` and `main_plot=TRUE`.

Notes:

- STAR points are taken from corrected HEPData tables using `alpha_Lambda = 0.732` where available.
- STAR 19.6 and 27 GeV BES-I points remain in the CSV, but the main plot uses the newer high-statistics BES-II values from Phys. Rev. C 108.
- STAR 200 GeV in the main plot uses the high-statistics Phys. Rev. C 98 point, not the older wide-uncertainty value in the Nature energy-scan table.
- ALICE main-plot points use the 15-50% wide centrality bin from the erratum-corrected HEPData record.
- HADES mid-central points are from the published Physics Letters B abstract; the detailed acceptance follows the paper.

The plotted vertical uncertainties are statistical and systematic uncertainties added in quadrature. Asymmetric systematic uncertainties are kept in the CSV and propagated asymmetrically.

`lambda_global_polarization_theory.csv` contains model curves digitized from the published figure images. HADES Fig. 3 uses `sqrt(sNN)-2mN` on the x-axis; the CSV stores the converted `sqrt_snn_GeV` values using `2mN = 1.876 GeV`.

`star_lambda_antilambda_global_polarization.csv` is a compact STAR-only extraction of the selected Lambda and anti-Lambda points used for the comparison figure.
