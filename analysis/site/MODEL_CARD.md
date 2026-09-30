# ReefDNA eDNA-first site-visit model

## Intended use

This is the current operational model. It predicts CPUE for the first
subsequent visit to each cull site following an eDNA sampling campaign.

The linkage is:

1. group same-reef eDNA events separated by no more than seven days into a
   sampling campaign;
2. identify the first subsequent visit to every cull site within 2,000 m and
   183 days;
3. select the nearest eDNA sampling site within that campaign;
4. resolve competing campaigns in favour of the most recent eligible campaign;
5. aggregate all cull dives at the chosen site and date by summing COTS counts
   and bottom time.

This produces 519 unique site visits across 332 sites, 78 eDNA campaigns and
36 reefs. No CPUE outcome is duplicated.

## Models

The BRMS formulation is:

    COTS count ~ % positive + distance + elapsed time
                 + offset(log(total bottom time))
                 + (1 | reef) + (1 | cull site) + (1 | eDNA campaign)

The concentration-free BRT uses the same three fixed predictors and a Poisson
count objective with log effort as its base margin. A design-aware experimental
BRT additionally includes indicators for `4x6` and `other_or_mixed`, with
`3x12` as the reference. The design-aware INLA screen compares the base model,
an additive sampling-design adjustment, and a sampling-design by percent-positive
interaction, while retaining reef, cull-site and campaign random intercepts.

The 519 response rows comprise 160 `3x12`, 126 `4x6`, and 233
`other_or_mixed` site visits. The latter category preserves older and transition
voyages rather than incorrectly assigning them to a principal protocol.

The existing production CPUE ensemble uses 59% BRMS and 41% BRT. Threshold
probabilities have separately calibrated weights. It should not be relabelled as
design-adjusted until its weights are re-estimated using the new BRT predictions.

## Validation

Five-fold validation holds out complete eDNA campaigns. The long production
BRMS model and all final validation fits have zero divergent transitions.

| Model | RMSE | MAE | Correlation |
|---|---:|---:|---:|
| Ensemble | 0.0471 | 0.0318 | 0.440 |
| BRMS | 0.0479 | 0.0329 | 0.436 |
| BRT | 0.0487 | 0.0322 | 0.363 |

At 0.04 CPUE, the ensemble has AUC 0.745, Brier score 0.180, recall 0.745 and
specificity 0.624 at its out-of-fold F1 cutoff. The cutoff is exploratory
because it was selected on the same out-of-fold predictions used to report F1.

In the design-aware screen, the INLA design-by-percent-positive formulation had
the lowest WAIC (3527.0), 7.5 below the additive-design model and 8.8 below the
base model. This in-sample improvement did not transfer to the BRT campaign-held-
out assessment: RMSE was 0.0489 with design versus 0.0487 without it. Sampling
design contributed about 0.2% of BRT gain. Treat the INLA interaction as a
hypothesis for validation, not as established predictive improvement.

## Prediction

    models <- readRDS(
      "analysis/site/output/reefDNA_site_visit_model_bundle.rds"
    )

    new_sites <- data.frame(
      perc_pos = c(50, 75, 95),
      distance_m = c(200, 500, 1500),
      lag_days = c(30, 90, 180),
      Reef = "Example Reef",
      cull_site_name = c("Site 1", "Site 2", "Site 3"),
      edna_campaign_id = "example_2026_09"
    )

    reefDNA::predict_reef_cpue(models, new_sites, method = "ensemble")

If reef, site, or campaign identifiers match training levels, their fitted
effects are used. New identifiers are assigned draws from the corresponding
hierarchical distribution. If bottom time is omitted, the median training
site-visit effort of 448 minutes is used for the observation interval and
threshold probabilities.

Outputs include expected CPUE, an expected-value interval, a future-observation
interval, probabilities of exceeding 0.02, 0.04 and 0.08 CPUE, and provisional
F1-tuned alert flags.

## Limitations

- Distance has an uncertain fixed effect after accounting for reef, site and
  campaign variability.
- Sampling design is strongly associated with voyage era and funding program;
  its effect is vulnerable to temporal and program confounding.
- The dataset contains only 36 reefs; predictions on new reefs carry wide
  hierarchical uncertainty.
- BRT intervals are empirical whole-reef bootstrap or out-of-fold residual
  intervals. The ensemble is pragmatic rather than a coherent joint posterior.
- Predictions beyond 2,000 m or 183 days are extrapolations outside the linkage
  design and should not be used operationally without additional validation.
