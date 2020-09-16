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
tinytest::expect_identical(
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
tinytest::expect_equal(
    MSstatsConvert:::.getOverlappingFeatures(data_with_additional_cols),
    c("A_B_3", "A_C_3", "A_D_3")
)
## Non-overlapped feature is not considered
tinytest::expect_equal(
    MSstatsConvert:::.filterOverlapped(data_with_additional_cols, summary_function = mean,
                                       MSstatsConvert:::.getOverlappingFeatures(data_with_additional_cols)),
    data_with_additional_cols[-(c(3:4, 19:20)), ]
)
