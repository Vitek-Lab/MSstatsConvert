# Data for testing fractions ----
## Simplest case:
data_for_fractions = data.table::data.table(
    ProteinName = "A",
    PeptideSequence = rep(rep(c("A", "B", "C", "D"), times = 1:4), each = 2),
    PrecursorCharge = 3,
    PSM = paste(rep(rep(c("A", "B", "C", "D"), times = 1:4), each = 2), 3, sep = "_"),
    Mixture = 1,
    TechRepMixture = 1,
    Run = as.character(rep(1:10, each = 2)),
    Channel = rep(c(1, 2), times = 2),
    Intensity = c(1, 1, 1, 1, 2, 2, 1, NA, 1, NA, 0, 2, 1, NA, 1, NA, 1, NA, -0.1, 1.1),
    Fraction = c(1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4)
)
# Tests ----
## No fractions - do nothing
expect_identical(
    MSstatsConvert:::.handleFractions(data_for_fractions[1:2, ]),
    data_for_fractions[1:2, ]
)
## We find the right features
data_with_additional_cols = data_for_fractions
data_with_additional_cols$techrun = paste(data_with_additional_cols$Mixture, 
                                          data_with_additional_cols$TechRepMixture, 
                                          sep = "_")
.combine = MSstatsConvert:::.combine
data_with_additional_cols$feature = do.call(".combine", list(data_with_additional_cols$ProteinName,
                                                             data_with_additional_cols$PeptideSequence,
                                                             data_with_additional_cols$PrecursorCharge))
data_with_additional_cols$id = paste(data_with_additional_cols$feature, 
                                     data_with_additional_cols$Run, sep = "_")
expect_equal(
    MSstatsConvert:::.getOverlappingFeatures(data_with_additional_cols),
    c("A_B_3", "A_C_3", "A_D_3")
)
## Non-overlapped feature is not considered
expect_equal(
    MSstatsConvert:::.filterOverlapped(data_with_additional_cols, summary_function = mean,
                                       MSstatsConvert:::.getOverlappingFeatures(data_with_additional_cols)),
    data_with_additional_cols[-(c(3:4, 19:20)), ]
)
## Overlap check
### Fails, when overlapped
expect_error(
    MSstatsConvert:::.checkOverlappedFeatures(
        data.table::data.table(
            feature = "A",
            Fraction = c(1, 2, 3),
            Intensity = c(NA, 1, 2)
        )
    )
)
### OK, when not overlapped
fractions_ok = data.table::data.table(
    feature = "A",
    Fraction = c(1, 2, 2),
    Intensity = c(NA, 1, 2)
)
expect_identical(
    MSstatsConvert:::.checkOverlappedFeatures(fractions_ok),
    fractions_ok
)
## 
fractionated = data.table::data.table(
    feature = rep(c("A", "B"), each = 6),
    Fraction = rep(rep(c(1, 2), each = 3), times = 2),
    Run = 1:12,
    IsotopeLabelType = "L",
    Intensity = c(NA, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2)
)
### More observations win
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated[feature == "A"])$Fraction),
    2
)
### Higher average intensity wins on tie
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated[feature == "B"])$Fraction),
    2
)
### For full data
expect_identical(
    MSstatsConvert:::.removeOverlappingFeatures(data.table::copy(fractionated)),
    fractionated[Fraction == 2]
)
### Non-tied fraction with high mean intensity should not be selected over tied fractions
fractionated_third = data.table::data.table(
    feature = rep("A", 9),
    Fraction = c(rep(1, 3), rep(2, 3), rep(3, 3)),
    Run = 1:9,
    IsotopeLabelType = "L",
    Intensity = c(1, 1, 1,    # Fraction 1: 3 obs, mean = 1 (ties for max n_obs)
                  2, 2, 2,    # Fraction 2: 3 obs, mean = 2 (ties for max n_obs)
                  10, NA, NA) # Fraction 3: 1 obs, mean = 10 (loses on n_obs but mean would win w/o fix)
)
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated_third)$Fraction),
    2
)
### L/H isotope label types: fraction selection uses only L obs, output keeps both L and H;
### NA-only features: fraction selection uses NA obs
# Feature "A" (L+H): F1 has 1 L obs (run 1) + 4 H obs (runs 1-4) = 4 unique runs total
#                    F2 has 3 L obs (runs 5-7) + 1 H obs (run 5)  = 3 unique runs total
#   → If counting all unique runs: F1 wins (4 > 3)
#   → If counting L unique runs only: F2 wins (3 > 1)
#   → Expected: F2 selected, confirming H obs do not influence fraction selection
#   → Both L and H rows from F2 are returned
# Feature "B" (NA only): F1 has 2 NA obs (runs 8-9), F2 has 4 NA obs (runs 10-13)
#   → NA obs: F1=2, F2=4 → F2 wins; all 4 F2 rows are kept
fractionated_lh = data.table::data.table(
    feature = c(rep("A", 9), rep("B", 6)),
    Fraction = c(rep(1, 5), rep(2, 4),
                 rep(1, 2), rep(2, 4)),
    Run = c(1, 1, 2, 3, 4,    # A F1: run 1 paired L+H, runs 2-4 H only
            5, 5, 6, 7,        # A F2: run 5 paired L+H, runs 6-7 L only
            8, 9, 10, 11, 12, 13),
    IsotopeLabelType = c("L", "H", "H", "H", "H",   # A F1: 1 L, 4 H
                         "L", "H", "L", "L",          # A F2: 3 L, 1 H
                         rep(NA_character_, 6)),
    Intensity = rep(1, 15)
)
### Feature A: F1 has more total obs but fewer L obs — F2 must win to confirm L-only logic
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated_lh[feature == "A"])$Fraction),
    2
)
### Both L and H rows are returned for the winning fraction
expect_true(
    all(c("L", "H") %in% MSstatsConvert:::.removeOverlappingFeatures(fractionated_lh[feature == "A"])$IsotopeLabelType)
)
### Feature B (NA only) selected via NA obs count
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated_lh[feature == "B"])$Fraction),
    2
)
### Full dataset: A keeps F2 (L+H rows), B keeps F2 (NA rows)
expect_identical(
    MSstatsConvert:::.removeOverlappingFeatures(data.table::copy(fractionated_lh)),
    fractionated_lh[(feature == "A" & Fraction == 2) | (feature == "B" & Fraction == 2)]
)
### H-only feature: falls back to H obs for fraction selection (not silently dropped)
# F1 has 1 H obs, F2 has 3 H obs → F2 wins
fractionated_h_only = data.table::data.table(
    feature = rep("A", 4),
    Fraction = c(1, 2, 2, 2),
    Run = 1:4,
    IsotopeLabelType = "H",
    Intensity = 1
)
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated_h_only)$Fraction),
    2
)
### H-only feature with tied obs count: higher mean H intensity breaks the tie
# F1: 2 obs, mean intensity 1; F2: 2 obs, mean intensity 3 → F2 wins
fractionated_h_tied = data.table::data.table(
    feature = rep("A", 4),
    Fraction = rep(c(1, 2), each = 2),
    Run = 1:4,
    IsotopeLabelType = "H",
    Intensity = c(1, 1, 3, 3)
)
expect_equal(
    unique(MSstatsConvert:::.removeOverlappingFeatures(fractionated_h_tied)$Fraction),
    2
)
fractionated_tmt = fractionated = data.table::data.table(
    feature = rep(c("A", "B"), each = 6),
    Fraction = rep(rep(c(1, 2), each = 3), times = 2),
    Run = rep(rep(c(1, 2), each = 3), times = 2),
    Intensity = c(NA, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2),
    Channel = rep(1:3, times = 4),
    Mixture = rep(c(1, 2), each = 6),
    TechRepMixture = 1
)
tmt_unoverlapped = MSstatsConvert:::.handleFractionsTMT(data.table::copy(fractionated_tmt))
expect_identical(
    tmt_unoverlapped[, .(Mixture, TechRepMixture, Channel, feature, Intensity)],
    fractionated_tmt[Fraction == 2, .(Mixture, TechRepMixture, Channel, feature, Intensity)]
)
