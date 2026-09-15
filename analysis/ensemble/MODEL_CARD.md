# ReefDNA operational CPUE model card

> Superseded for operational use by the eDNA-first site-visit model in
> analysis/site/MODEL_CARD.md. This object is retained as a cull-first
> sensitivity analysis and for reproducibility.

## Purpose

The model bundle predicts crown-of-thorns starfish CPUE from eDNA percent
positive, distance from the eDNA site, and elapsed time since eDNA sampling.
Concentration is deliberately excluded from the operational formulation.

## Model bundle

Build the local bundle with:

```r
source("analysis/ensemble/build_model_bundle.R")
```

This creates the git-ignored file
`analysis/ensemble/output/reefDNA_operational_model_bundle.rds`. It contains:

- the final negative-binomial BRMS model;
- the final Poisson boosted-tree regression and three threshold classifiers;
- predictor scaling rules;
- whole-reef cross-validation ensemble weights; and
- F1-tuned alert cutoffs.

## Required prediction inputs

Supply a data frame with one row per prediction and these columns:

| Input | Meaning | Accepted aliases |
|---|---|---|
| `perc_pos` | eDNA-positive samples, 0 to 100 percent | `edna_percent_positive`, `edna_pct` |
| `distance_m` | distance from the eDNA site in metres | `edna_distance_m` |
| `lag_days` | days since eDNA sampling | `edna_lag_days` |
| `bottom_time` | survey effort in minutes; defaults to 216 | `survey_minutes` |
| `Reef` | optional reef name | none |

Example:

```r
models <- readRDS(
  "analysis/ensemble/output/reefDNA_operational_model_bundle.rds"
)

new_edna <- data.frame(
  perc_pos = c(40, 75, 95),
  distance_m = c(200, 500, 1500),
  lag_days = c(30, 90, 180),
  bottom_time = 216
)

reefDNA::predict_reef_cpue(models, new_edna, method = "ensemble")
```

The result gives expected CPUE, uncertainty and future-observation intervals,
probabilities of exceeding 0.02, 0.04 and 0.08 CPUE, and optional alert flags.

## Validation summary

All reported validation predictions held out complete reefs in five folds.

| Model | RMSE | MAE | Correlation |
|---|---:|---:|---:|
| Ensemble | 0.0597 | 0.0393 | 0.298 |
| BRT | 0.0603 | 0.0397 | 0.302 |
| BRMS | 0.0618 | 0.0404 | 0.177 |
| GLMM | 0.0622 | 0.0401 | 0.181 |

The CPUE ensemble uses 34% BRMS and 66% BRT. Threshold probabilities use
separately tuned weights because their calibration differs by threshold.

## Interpretation and limitations

- The earlier low BRMS curve was mainly a prediction-scale error: setting the
  reef random effect to zero estimates a median or typical reef. New-reef
  predictions now integrate between-reef variation.
- The concentration-by-distance interaction shifts the modelled 0.04 crossing
  later (about 69% to 75% positive), but was not the main reason the earlier
  curve failed to cross 0.04.
- Only 36 reefs contributed to training. Whole-reef discrimination is modest,
  so outputs are decision support rather than a stand-alone deployment rule.
- F1-tuned alert cutoffs strongly favour sensitivity. At 0.04, the ensemble
  cutoff is 0.07 probability and has approximately 99.9% recall but 1.2%
  specificity on the same out-of-fold predictions used to choose it. Use the
  calibrated probability columns by default; deploy alerts only after agreeing
  the relative costs of false negatives and false positives.
- BRT intervals use out-of-fold residual quantiles. Ensemble intervals combine
  those empirical limits with BRMS posterior intervals and should be treated as
  pragmatic operational intervals, not a single coherent posterior.

The complete analysis is in `01_operational_model_assessment.html` in this
directory.
