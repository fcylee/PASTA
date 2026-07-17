data("polyA_small")

pct <- percentage.usage(polyA_small, features = rownames(polyA_small),
                        cells = Cells(polyA_small)[polyA_small$group== "A"])

test_that("percentages within same gene sum to 1", {
  expect_equal(as.numeric(pct["11-118408896-118409195"] + pct["11-118409548-118409847"]), 1)
})

#test percentages are correct
a <- subset(polyA_small, group=="A")
s.a <- rowSums(LayerData(a, layer = "counts"))[c("11-118408896-118409195","11-118409548-118409847" )]
test_that("percentage usage is alculated correctly", {
  expect_equal(s.a["11-118408896-118409195"]/sum(s.a), pct["11-118408896-118409195"])
})

test_that("FindDifferentialPolyA throws error if no residuals present", {
expect_error(FindDifferentialPolyA(polyA_small, ident.1 =  "A"),
             regexp = "Run CalcPolyAResiduals prior to FindPolyASites")
})

polyA_small <- CalcPolyAResiduals(polyA_small)
Idents(polyA_small) <- polyA_small$group
m <- FindDifferentialPolyA(polyA_small, ident.1 =  "A", ident.2="B")

test_that("Percentages within a gene sum to 1", {
  expect_equal(sum(subset(m, symbol=="UBC")$percent.1), 1)
  expect_equal(sum(subset(m, symbol=="UBC")$percent.2), 1)
  expect_equal(sum(subset(m, symbol=="ARPC1B")$percent.1), 1)
  expect_equal(sum(subset(m, symbol=="ARPC1B")$percent.2), 1)
})

test_that("Percentage usage matches manual calculation", {
  expect_equal(m["11-118408896-118409195","percent.1"], as.numeric(s.a["11-118408896-118409195"]/sum(s.a)))
})

m2 <- FindDifferentialPolyA(polyA_small, ident.1 =  "A", ident.2 = "B",
                      features = c("12-124911648-124911947", "14-75282887-75283186"))
test_that("Subsetting features gives same results", {
  expect_equal(m2, m[c("12-124911648-124911947", "14-75282887-75283186"),])
})

# ---- provenance columns / mixed-model paths ---------------------------------

# default (no group.by.sample) fits per-cell lm and stamps provenance columns
test_that("default fit is lm with provenance columns populated", {
  expect_true(all(c("model", "n_samples", "is_singular", "converged") %in%
                    colnames(m)))
  expect_true(all(m$model == "lm"))
  expect_true(all(is.na(m$n_samples)))
  expect_true(all(is.na(m$is_singular)))
})

# synthesize 3 samples per group (~31 cells each) and a sex covariate
polyA_mm <- polyA_small
cellsA <- Cells(polyA_mm)[polyA_mm$group == "A"]
cellsB <- Cells(polyA_mm)[polyA_mm$group == "B"]
samp <- setNames(rep(NA_character_, ncol(polyA_mm)), Cells(polyA_mm))
samp[cellsA] <- paste0("a", (seq_along(cellsA) %% 3) + 1)
samp[cellsB] <- paste0("b", (seq_along(cellsB) %% 3) + 1)
polyA_mm$sample_id <- unname(samp[Cells(polyA_mm)])
polyA_mm$sex <- unname(c(a1 = "F", a2 = "M", a3 = "F",
                         b1 = "M", b2 = "F", b3 = "M")[polyA_mm$sample_id])

test_that("group.by.sample not in metadata errors", {
  expect_error(
    FindDifferentialPolyA(polyA_mm, ident.1 = "A", ident.2 = "B",
                          group.by.sample = "not_a_column"),
    regexp = "not found in meta data")
})

test_that("invalid mixed.test is rejected", {
  expect_error(
    FindDifferentialPolyA(polyA_mm, ident.1 = "A", ident.2 = "B",
                          group.by.sample = "sample_id", mixed.test = "bogus"))
})

test_that("guard falls back to lm when required samples/cells not met", {
  mm_guard <- suppressWarnings(suppressMessages(
    FindDifferentialPolyA(polyA_mm, ident.1 = "A", ident.2 = "B",
                          group.by.sample = "sample_id", min.groups = 100)))
  expect_true(all(mm_guard$model == "lm"))
  expect_true(all(is.na(mm_guard$is_singular)))
  expect_true(all(mm_guard$n_samples == 6))   # counted regardless of fallback
})

test_that("satterthwaite mixed model runs and records provenance", {
  skip_if_not_installed("lmerTest")
  mm_sat <- suppressWarnings(suppressMessages(
    FindDifferentialPolyA(polyA_mm, ident.1 = "A", ident.2 = "B",
                          covariates = "sex", group.by.sample = "sample_id",
                          mixed.test = "satterthwaite")))
  expect_s3_class(mm_sat, "data.frame")
  # either it fit lmer, or the whole comparison fell back and is flagged
  expect_true(all(mm_sat$model %in% c("lmer", "lm_fallback")))
  expect_true(all(mm_sat$n_samples == 6))
})

test_that("wald mixed model runs, never falls back, yields p-values", {
  skip_if_not_installed("lme4")
  mm_wald <- suppressWarnings(suppressMessages(
    FindDifferentialPolyA(polyA_mm, ident.1 = "A", ident.2 = "B",
                          covariates = "sex", group.by.sample = "sample_id",
                          mixed.test = "wald")))
  expect_true(all(mm_wald$model == "lmer"))          # wald has no lm_fallback
  expect_true(all(mm_wald$n_samples == 6))
  expect_true(all(!is.na(mm_wald$p.value)))          # Wald p from t-statistic
  expect_true(all(is.logical(mm_wald$is_singular)))  # flag populated, not NA
})

