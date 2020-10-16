# Example data ----
## Without channels
to_filter = data.table::data.table(
    PeptideSequence = sort(rep(letters[1:20], times = 3)),
    PrecursorCharge = rep(rep(c(2, 3), each = 15), times = 2),
    Run = rep(rep(1:3, times = 5), times = 2),
    Intensity = seq(0.1, 10, length.out = 60)
)
## With channels
to_filter_channel = data.table::data.table(
    ProteinName = "A",
    PeptideSequence = sort(rep(rep(letters[1:10], times = 3), each = 5)),
    PrecursorCharge = 3,
    Run = rep(rep(1:3, each = 5), 10),
    Channel = rep(1:5, times = 30),
    Fraction = 1,
    Intensity = seq(0.1, 10, length.out = 150)
)
to_filter_channel$Intensity[1:6] = NA
to_filter_channel$PSM = paste(to_filter_channel$PeptideSequence,
                              to_filter_channel$PrecursorCharge,
                              sep = "_")
# Remove feature with few measurements ----
# Filtering for TMT
## Filter few works for TMT data - removes PSMs with less than 3 observations
## across channels
expect_equal(
    MSstatsConvert:::.filterFewMeasurements(
        data.table::copy(to_filter_channel[1:15]), 0, TRUE, c("PSM", "Run")),
    to_filter_channel[(PeptideSequence == "a") & Run %in% c(2, 3), ]
)
## Nothing is missing
expect_equal(
    nrow(MSstatsConvert:::.filterFewMeasurements(
        data.table::copy(to_filter)[, Intensity := 0.1], 0, TRUE, c("PeptideSequence", "PrecursorCharge"))),
    nrow(to_filter)
)
## With a threshold:
expect_equal(
    MSstatsConvert:::.filterFewMeasurements(
        data.table::copy(to_filter), 1, TRUE, c("PeptideSequence", "PrecursorCharge")),
    to_filter[!(PeptideSequence %in% c("a", "b"))]
)
## With missing values
to_filter$Intensity = ifelse(to_filter$PeptideSequence %in% c("a", "b", "c") &
                                 to_filter$Run %in% 1:2, NA, to_filter$Intensity)
to_filter$Intensity[3] = NA
## Remove only features that have all missing values
expect_equal(
    MSstatsConvert:::.filterFewMeasurements(
        data.table::copy(to_filter), 0, FALSE, c("PeptideSequence", "PrecursorCharge")),
    to_filter[PeptideSequence != "a", ]
)
## Remove all features with less than three non-missing values
expect_equal(
    MSstatsConvert:::.filterFewMeasurements(
        data.table::copy(to_filter), 0, TRUE, c("PeptideSequence", "PrecursorCharge")),
    to_filter[!(PeptideSequence %in% c("a", "b", "c"))]
)
# Aggregating / summarizing features ----
to_aggregate = data.table::data.table(PeptideSequence = "A", 
                                      Run = rep(1:3, each = 2),
                                      Intensity = 1:6)
to_aggregate$Intensity[5] = NA
expect_equal(
    MSstatsConvert:::.summarizeMultipleMeasurements(to_aggregate, max, c("PeptideSequence", "Run")),
    to_aggregate[c(2, 4, 6), ]
)
## Zeros are coded correctly
### 
with_zeros = data.table::data.table(
    ProteinName = "A",
    PeptideSequence = "a",
    FragmentIon = rep(c("precursor", "precursor [M+1]", "precursor [M+2]"),
                      times = 4),
    Run = rep(1:4, each = 3),
    Intensity = c(0, 2, 3, NA, NA, NA, 0, 0, NA, 0, 0, 0),
    isZero = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE,
               TRUE, TRUE, TRUE)
)
correct_result = data.table::data.table(
    PeptideSequence = "a",
    FragmentIon = NA,
    ProductCharge = NA,
    Run = 1:4,
    Intensity = c(5, 0, 0, 0),
    isZero = c(FALSE, FALSE, TRUE, TRUE),
    ProteinName = "A"
)
MSstatsConvert:::.checkDDA(with_zeros)
expect_equal(
    MSstatsConvert:::.summarizeMultipleMeasurements(with_zeros, sum, c("PeptideSequence", "Run")),
    correct_result[, -(2:3), with = FALSE]
)
expect_equal(
    MSstatsConvert:::.handleIsotopicPeaks(with_zeros, TRUE),
    correct_result[, c(7, 1, 4, 5, 6, 2, 3), with = FALSE]
)
expect_equal(
    MSstatsConvert:::.handleIsotopicPeaks(with_zeros, FALSE),
    with_zeros
)
## TMT data
to_aggregate_channel = data.table::copy(to_filter_channel)
to_aggregate_channel$Run = rep(1:6, each = 25)
to_aggregate_channel$PSM = paste(to_aggregate_channel$PeptideSequence,
                                 to_aggregate_channel$PrecursorCharge,
                                 rep(1:30, each = 5), sep = "_")
to_aggregate_channel$Intensity[1:15] = 1
to_aggregate_channel$Feature = MSstatsConvert:::.combine(to_aggregate_channel$PeptideSequence,
                                                         to_aggregate_channel$PrecursorCharge)
to_aggregate_channel[, n_psms := data.table::uniqueN(PSM), 
                     by = c("PeptideSequence", "PrecursorCharge")]
## PSM summarization works as expected
### Case 1: there is a score
to_aggregate_channel$Score = 1
to_aggregate_channel$Score[16:30] = rep(1:3, each = 5)
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[16:30, ], max),
    "b_3_6"
)
### Case 3: there is isolation interference
to_aggregate_channel$IsolationInterference = 1
to_aggregate_channel$IsolationInterference[31:45] = rep(1:3, each = 5)
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[31:45, ], max),
    "c_3_7"
)
### Case 4: there is IonsScore
to_aggregate_channel$IonsScore = 1
to_aggregate_channel$IonsScore[46:60] = rep(1:3, each = 5)
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[46:60, ], max),
    "d_3_12"
)
### Case 5: number of nonmissing values
to_aggregate_channel$Intensity[c(61:63, 66:67, 71)] = NA 
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[61:75, ], max),
    "e_3_15"
)
### Case 6: maximal intensity
to_aggregate_channel$Intensity[76:90] = 1:15
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[76:90, ], max),
    "f_3_18"
)
### Case 7: cannot decide - average intensities
to_aggregate_channel2 = data.table::copy(to_aggregate_channel[PeptideSequence == "b"][1:10])
to_aggregate_channel2$Intensity = 1
to_aggregate_channel2[, IsolationInterference := NULL]
to_aggregate_channel2[, IonsScore := NULL]
to_aggregate_channel2[, Score := NULL]
expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel2, max),
    NA
)
to_compare_agg_psm = to_aggregate_channel2[1:10, list(Intensity = mean(Intensity, na.rm = TRUE),
                                                      PSM = "b_3"),
                                           by = c("ProteinName", "PeptideSequence", "PrecursorCharge",
                                                  "Run", "Channel", "Fraction", "Feature")]
# to_compare_agg_psm$PSM = to_compare_agg_psm$Feature
expect_equal(
    MSstatsConvert:::.aggregatePSMstoPeptideIons(
        data.table::copy(to_aggregate_channel2[1:10, ]),
        c("PeptideSequence", "PrecursorCharge")),
    to_compare_agg_psm
)
## Aggregation works as expected
aggregated = MSstatsConvert:::.aggregatePSMstoPeptideIons(to_aggregate_channel, c("PeptideSequence", "PrecursorCharge"), summary_function = sum)
expect_equal(nrow(aggregated), 5 * (data.table::uniqueN(to_aggregate_channel[, .(PeptideSequence, PrecursorCharge, Run)])))
# Removing proteins with just a single feature ----
to_filter_channel2 = to_filter_channel[, !(colnames(to_filter_channel) %in% c("feature", "feature_count")), 
                                       with = FALSE]
dataset_tmt2 = rbind(to_filter_channel2, 
                     data.table::data.table(ProteinName = "B", PeptideSequence = "z",
                                            PSM = "z_3", PrecursorCharge = 3, Run = rep(1:3, each = 5),
                                            Fraction = 1, Channel = rep(1:5, times = 3), Intensity = 1))
## Do nothing
expect_equal(
    MSstatsConvert:::.handleSingleFeaturePerProtein(to_filter_channel2, FALSE),
    to_filter_channel2
)
## Filter when needed
expect_equal(
    MSstatsConvert:::.handleSingleFeaturePerProtein(to_filter_channel2, TRUE),
    data.table::setkey(to_filter_channel2[!(ProteinName == "B"), ], NULL)
)

# Utility function ----
expect_equal(MSstatsConvert:::.combine(c("A", "B"), c("A", "B")),
             c("A_A", "B_B"))
# Top-level function works
clean_lf = MSstatsConvert:::.cleanByFeature(
    to_filter, 
    c("PeptideSequence", "PrecursorCharge"),
    list(remove_features_with_few_measurements = TRUE,
         summarize_multiple_psms = max))
clean_tmt = MSstatsConvert:::.cleanByFeature(
    to_filter_channel, 
    c("PeptideSequence", "PrecursorCharge"),
    list(remove_features_with_few_measurements = TRUE,
         summarize_multiple_psms = max))
expect_equal(class(clean_lf)[1], "data.table")
expect_equal(class(clean_tmt)[1], "data.table")

