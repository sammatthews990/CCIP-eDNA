# Spatio-temporal INLA plan for eDNA-informed COTS CPUE

## Objective

Fit a Bayesian spatio-temporal model that uses the raw eDNA replicate observations (binary detection and concentration), their sampling locations and dates, and spatial/temporal separation from cull observations to estimate expected COTS CPUE and its uncertainty. The required products are:

- site-level posterior CPUE at 0, 200, 500, 1000, and 2000 m; and at 3. 6 and 12 months
- reef-level, area-weighted posterior CPUE;
- 50%, 80%, and 95% credible intervals;
- posterior probabilities that CPUE exceeds 0.02 and 0.04; and
- maps of posterior mean, uncertainty, and management-threshold exceedance probability.

The final model should retain raw eDNA replicates in their own observation likelihood. It should not copy the same cull response onto every eDNA replicate, because that would treat repeated laboratory or field samples as independent CPUE observations and produce over-confident inference.

## Implementation status (7 September 2026)

Stage A is now implemented and fitted in `analysis/inla/01_spatiotemporal_benchmark.qmd` using reusable functions in `reefDNA/R/inla_data.R`, `reefDNA/R/inla_mesh.R`, and `reefDNA/R/inla_models.R`.

- The audited linkage retains 7,612 raw eDNA rows and 67,826 admissible raw links, then collapses them to 2,117 unique cull-dive responses across 36 reefs before modelling.
- Static cull-site polygons are linked by normalized name first and a documented same-reef nearest-polygon fallback within 2,000 m. This leaves three candidate cull rows unmatched in the model date window.
- The fitted negative-binomial benchmark uses effort as an offset, reef heterogeneity, an annual AR(1) Matérn field, and separate eDNA, distance, and lag covariates.
- It produces 785 predictions for 157 eDNA sites at 0, 200, 500, 1,000, and 2,000 m, with 80% and 95% credible intervals and threshold exceedance probabilities.
- A lighter operational screen is implemented in `analysis/inla/02_distance_time_formulation_screen.qmd`. It pairs each cull with its most recent eligible prior eDNA event, then compares all 16 additive subsets of percent positive, mean concentration, distance from the sample, and days since the sample, plus seven hierarchical pairwise-interaction formulations. It excludes the SPDE and calendar year while retaining a reef iid intercept.
- The best additive formulation is percent positive + distance + lag. The all-pairwise interaction model has the lowest WAIC, but its scenario surface reverses the distance effect at long lags; retain both for grouped cross-validation rather than accepting the interaction model on WAIC alone.
- A `brms` negative-binomial validation of this shortlist, using grouped folds by reef or eDNA event and posterior predictive checks, is now the next practical milestone. Stage B, the shared-field joint Bernoulli/detected-concentration/count model, remains a later scientific extension rather than the default operational model.

## Data-readiness snapshot (7 September 2026)

- `data/eDNA data_ALL_20260528.xlsx`: 7,612 rows, all with date and coordinates; 3,952 detections; 3,947 positive concentrations; 720 distinct reef-site-year combinations; observations span 6 May 2019 to 8 April 2026.
- `data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx`: 49,130 rows with valid date, coordinates, bottom time, and positive effort; observations span 22 February 2002 to 29 May 2026. There are 26,680 usable dives within the eDNA date range.
- `data/Eotr_CotsCullSites_2025_11_19_1_58_PM.gpkg`: 15,884 static cull-site multipolygons in EPSG:4326. It has site names but no survey dates or effort.
- `data/rrap_canonical_2025-03-20-T15-18-17.gpkg`: 3,806 reef polygons in EPSG:7844.
- The initial name-only scan linked 88.1% of unique raw cull-site names but only 59.0% of valid cull-dive rows. The implemented full reconciliation uses normalized name matching plus a same-reef nearest-polygon fallback within 2,000 m; among 27,104 candidate rows it links 27,042 by name, 59 by fallback, and leaves three unmatched.
- All distances, buffers, and mesh operations should use EPSG:3112. Convert to EPSG:4326 only for web tiles and longitude/latitude presentation.

## Required modelling table

Create three linked tables rather than one duplicated join:

1. `edna_obs`: one row per raw eDNA sample, with stable sample ID, reef, site, date, year/quarter index, coordinates, collection organisation, detection, concentration, PCR values where useful, and sampling-design metadata.
2. `cull_obs`: one row per cull dive, with dive ID, reef, cull-site name, date, coordinates, bottom time, total COTS count, cohort counts, depth, vessel/source, and polygon link status.
3. `space_time_links`: audited links between eDNA observations and candidate cull observations/sites, including Euclidean distance in metres, signed lag in days, same-reef indicator, link method, and whether the pair falls inside the prespecified temporal and 2000 m spatial windows.

Use the raw cull-dive point as the response location. Use the polygon distance as the main operational distance where a cull-site polygon is linked; otherwise use distance to the recorded dive point and flag the fallback. Keep signed temporal lag separate from spatial distance. The primary actionability window should be eDNA sampled 0–183 days before the cull dive, with 0–91, 92–183, and 184–365 day sensitivity analyses.

Before modelling, produce a reconciliation report containing match rates by year, reef, vessel/source, and cull-site naming pattern; distributions of distance and lag; duplicated IDs; impossible coordinates; zero/negative effort; and the number of reefs and time periods contributing to each model.

## Model architecture

Use a staged approach so that data-linkage and spatial assumptions can be validated before fitting the full joint model.

### Stage A: dive-level benchmark

Fit cull count at dive location and date using a negative-binomial likelihood with log bottom time as an offset:

`COTS_count_j ~ NegBin(mu_j, theta)`

`log(mu_j) = log(bottom_time_j) + beta0 + beta_edna * eDNA_summary_j + f_distance(d_j) + f_lag(lag_j) + reef_j + vessel_j + x(s_j, t_j)`

For this benchmark only, summarize raw eDNA observations within a reef-site sampling event using detection proportion and a robust concentration summary. This model provides a stable baseline and exposes linkage or mesh problems, but is not the final raw-replicate model.

### Stage B: joint raw-replicate model (primary model)

Represent local COTS abundance by a shared latent spatio-temporal field `x(s,t)`. Use three likelihoods tied to that field:

- Cull counts: negative binomial with `log(bottom_time)` offset.
- eDNA detection: Bernoulli with logit link.
- Positive eDNA concentration: Gaussian on `log1p(concentration)` initially; compare Gamma/lognormal-style alternatives if residual diagnostics require them.

Treat non-detect concentration as part of a hurdle model: the detection likelihood models zero/non-detect probability, while the conditional concentration likelihood uses detected samples only. Include a sample-event or site-date iid effect so replicate rows inform measurement error without pretending to be new spatial abundance observations.

The shared-field coefficients estimate how detection and concentration scale with latent COTS abundance. Implement the likelihoods through separate INLA stacks and a shared/copy latent component. Do not interpret the association as causal unless sampling design and current/transport covariates support that claim.

## Distance and time effects

The phrase “distance to cull site” can refer to two different quantities and they should not be silently mixed:

- Planned-site distance: distance from an eDNA location or prediction point to the nearest static GPKG cull polygon. This is available for prospective operational predictions.
- Observed-dive distance: distance to the nearest dated cull dive satisfying the reef and lag rules. This is useful for retrospective validation but may not exist at deployment time.

Use planned-site distance in the primary operational model and observed-dive distance as a sensitivity analysis. Start with `log1p(distance_m / 200)` for a stable, interpretable effect. Then compare a first-order random walk over prespecified distance bins or a low-complexity spline. Any nonlinear curve must be constrained by the observed support; do not extrapolate beyond 2000 m for the requested products.

Use a signed lag term rather than only broad windows. Start with `log1p(lag_days)` or a low-complexity RW1 term, then test a distance-by-lag interaction only after the additive model is stable. Add sampling year/season as fixed or cyclic effects if coverage supports them.

## Spatial and temporal latent structure

1. Build a two-dimensional mesh in EPSG:3112 over the sampled GBR domain, with an outer extension to reduce boundary effects. Prevent implausible short connections across unsampled land/open-water gaps using a non-convex or barrier-domain sensitivity model.
2. Use a Matérn SPDE created with `inla.spde2.pcmatern()` and explicit penalised-complexity priors for spatial range and marginal standard deviation.
3. Project cull, eDNA, and prediction locations with `inla.spde.make.A()`.
4. Couple spatial fields over time using an AR(1) group model. Begin with year or survey-quarter groups, depending on an empirical coverage table; a sparse monthly field should not be the default.
5. Retain a reef iid intercept to capture persistent reef-level heterogeneity not represented by the continuous field. Test site-event and vessel/source effects where identifiable.

Initial priors should be elicited on interpretable scales and then stress-tested. A defensible starting sensitivity grid is `P(range < 10, 25, or 50 km) = 0.05` and `P(spatial SD > 0.5, 1, or 2 on the link scale) = 0.05`. Set fixed-effect priors after centring/scaling continuous predictors. Document priors in the report and show whether management outputs materially change under the sensitivity grid.

## Candidate model sequence

Fit the following sequence and stop adding complexity when held-out predictive performance and calibration no longer improve:

1. Negative-binomial count plus effort offset and reef iid effect.
2. Add eDNA event summaries, distance, and lag.
3. Add spatial SPDE.
4. Add grouped AR(1) spatio-temporal field.
5. Replace event summaries with the joint Bernoulli plus conditional-concentration measurement model.
6. Add nonlinear distance/lag effects and optional covariates only if supported.

Useful optional covariates include depth, collection organisation, sampling design, season, and cull source/vessel. Current speed/direction, tide, and time since culling would be scientifically valuable for eDNA transport, but they should be listed as unavailable dependencies until data are supplied.

## Prediction products

Create a prediction stack that does not contain observed CPUE responses. For each requested date or scenario:

- predict on eDNA sites at distance scenarios 0, 200, 500, 1000, and 2000 m;
- predict across reef polygons on a resolution justified by the mesh;
- convert expected count back to CPUE at a fixed reference effort;
- retain posterior mean, median, SD, and 50%, 80%, and 95% credible intervals;
- calculate `Pr(CPUE > 0.02)` and `Pr(CPUE > 0.04)`; and
- aggregate grid cells to reefs with area weights, reporting both the posterior of reef-average CPUE and within-reef heterogeneity.

For operational communication, classify posterior probabilities rather than posterior means alone. The 0.02–0.04 band should be retained as a worthwhile/near-miss tier, consistent with the revised confusion matrix.

## Validation and acceptance criteria

Random row-wise cross-validation is not acceptable because replicates, sites, reefs, and nearby dates are correlated. Use:

- leave-one-reef-out or spatial-block cross-validation for geographic transfer;
- leave-one-year/campaign-out validation for temporal transfer;
- a combined space-time holdout as the most realistic operational test; and
- Eyrie Reef as a diagnostic case study, not as the only validation reef.

Compare against the existing negative-binomial GLMM and a non-spatial INLA model. Report held-out log score, RMSE/MAE for CPUE, interval coverage and width, calibration of `Pr(CPUE > threshold)`, sensitivity/precision/F1 at 0.02 and 0.04, and the three-tier near-miss table. Within-sample diagnostics should include CPO/PIT, WAIC, posterior predictive zero and upper-tail checks, dispersion, and maps of residual spatial autocorrelation.

Acceptance gates:

1. Data linkage: IDs are unique; all distances use EPSG:3112; polygon/direct-point fallback rates are reported; no cull response is multiplied by eDNA replicate count.
2. Mesh: no obvious cross-domain shortcuts; results are stable to reasonable mesh refinement.
3. Calibration: nominal 80% and 95% intervals have acceptable held-out coverage, and exceedance probabilities are reliability-calibrated.
4. Robustness: conclusions at 0.02 and 0.04 do not reverse under reasonable prior, time-bin, distance-definition, and likelihood choices.
5. Reproducibility: a clean session can rebuild derived data, fit the selected model, and regenerate tables/maps with fixed seeds and recorded INLA version.

## Implementation work packages

### 1. Data audit and linkage

- Create stable eDNA observation, sample-event, and cull-dive IDs.
- Normalize reef and cull-site names; explicitly review high-frequency unmatched patterns.
- Link cull rows to polygons by normalized name, then by same-reef spatial nearest neighbour under a documented tolerance.
- Generate coverage, distance, lag, and missingness reports.
- Freeze a versioned derived model dataset and data dictionary.

Deliverable: `analysis/inla/01_data_linkage.qmd` plus derived parquet/RDS tables and an audit CSV.

### 2. Mesh and benchmark model

- Build and plot the EPSG:3112 mesh, projector matrices, and time groups.
- Fit models 1–4 in the candidate sequence.
- Select time resolution and prior scale using coverage and block validation.

Deliverable: `analysis/inla/02_mesh_and_benchmark.qmd` and serialized mesh/model configuration.

### 3. Joint likelihood

- Build aligned response matrices/stacks for counts, detection, and positive concentration.
- Share the latent field with estimable scaling parameters.
- Add replicate/site-event effects and fit the hurdle measurement process.
- Check identifiability against simpler models and simulated data.

Deliverable: `analysis/inla/03_joint_model.qmd` and a model-comparison table.

### 4. Validation and sensitivity

- Run spatial, temporal, and combined block holdouts.
- Repeat the fit across prior, mesh, lag-window, distance-definition, and response-family choices.
- Choose the final model using predictive calibration and parsimony.

Deliverable: `analysis/inla/04_validation.qmd` and a locked model decision record.

### 5. Prediction and reporting

- Generate site-distance scenario curves and reef prediction surfaces.
- Produce reef summaries, credible intervals, exceedance probabilities, and three-tier management classifications.
- Add an Eyrie Reef map showing posterior predictions alongside the 200/500/1000/2000 m rings.

Deliverable: `analysis/inla/05_predictions.qmd`, geospatial outputs, and CSV tables for downstream use.

## Proposed package structure

- `reefDNA/R/inla_data.R`: validation, ID construction, name reconciliation, and space-time links.
- `reefDNA/R/inla_mesh.R`: domain, mesh, SPDE, projector, and prior helpers.
- `reefDNA/R/inla_models.R`: benchmark and joint model builders.
- `reefDNA/R/inla_predict.R`: site scenarios, reef integration, exceedance probabilities, and tidy posterior outputs.
- `reefDNA/tests/testthat/test-inla-data.R`: no duplicated response, CRS, lag sign, link fallback, and effort checks.
- `reefDNA/tests/testthat/test-inla-predict.R`: deterministic dimensions, thresholds, and area weights on a small simulated fit.
- `analysis/inla/config.yml`: temporal window, time resolution, mesh settings, priors, seeds, and model version.

Keep expensive fits out of routine package checks. Test data assembly and prediction functions with a small simulated mesh, and cache production fits with explicit hashes of input files and configuration.

## Decisions required before the production fit

These choices should be confirmed after the linkage audit rather than assumed now:

1. Is the operational predictor distance to a planned static cull polygon, a dated observed dive, or both?
2. Must eDNA always precede the CPUE observation, or should contemporaneous plus/minus windows also be modelled?
3. Is the production temporal resolution year, quarter, or campaign?
4. Should the primary reef prediction represent reef-wide average CPUE, expected CPUE at planned cull sites, or both?
5. Which collection/PCR metadata are scientifically meaningful batch effects?

## Recommended first executable milestone

Implement work package 1 and the Stage A benchmark first. This produces an auditable answer quickly, exposes spatial-linkage defects, and supplies empirical scales for the mesh and priors. Proceed to the joint raw-replicate model only after confirming that the benchmark can recover held-out reef/year CPUE without leakage.

## Updated operational modelling sequence

The completed screen changes the immediate priority. Use the broad SPDE/AR(1) model as a sensitivity benchmark, not the default operational formulation. The next `brms` analysis should fit the same one-event-per-cull table and compare this fixed shortlist:

1. percent positive + distance + lag (best additive INLA model);
2. all four additive main effects;
3. all main effects + concentration-by-distance;
4. all main effects + all pairwise interactions (INLA WAIC winner); and
5. the reef-intercept null model.

Use `neg_binomial_2(link = "log")`, `offset(log(bottom_time))`, standardized transformed predictors, and `(1 | Reef)` in every model. Validate with folds grouped by reef and, separately, by eDNA event/campaign. Compare held-out log predictive density, CPUE MAE/RMSE, interval coverage, and Brier/reliability scores for 0.02 and 0.04. Prefer the additive formulation unless the interaction model improves grouped prediction and preserves scientifically plausible distance-by-time surfaces. Add vessel/source effects as a sensitivity model; do not restore a year field or SPDE unless residual diagnostics show remaining structured dependence.

## R-INLA references

- SPDE overview and stack workflow: https://inla.r-inla-download.org/r-inla.org/doc/vignettes/SPDEhowto.pdf
- Matérn SPDE with PC priors: https://www.r-inla.org/learnmore/docs/reference/inla.spde2.pcmatern.html
- Observation/prediction matrices and temporal groups: https://www.r-inla.org/learnmore/docs/reference/inla.spde.make.A.html
- Group models including AR(1): https://www.r-inla.org/learnmore/docs/reference/control.group.html
