test_that("voyage sampling design distinguishes 3x12, 4x6 and mixed designs", {
  make_rows <- function(voyage, reef, n_sites, n_reps) {
    expand.grid(
      Voyage = voyage,
      ReefName = reef,
      Site_name = paste0("site_", seq_len(n_sites)),
      rep = seq_len(n_reps),
      stringsAsFactors = FALSE
    )
  }
  edna <- dplyr::bind_rows(
    make_rows("three", "reef_a", 3, 12),
    make_rows("four", "reef_b", 4, 6),
    make_rows("mixed", "reef_c", 3, 12),
    make_rows("mixed", "reef_d", 4, 6)
  )

  result <- summarise_edna_voyage_design(edna)
  expect_equal(result$sampling_design[match("three", result$Voyage)], "3x12")
  expect_equal(result$sampling_design[match("four", result$Voyage)], "4x6")
  expect_equal(result$sampling_design[match("mixed", result$Voyage)], "other_or_mixed")
})
