# Example data ----
## Without channels
to_filter = data.table::data.table(
    PeptideSequence = sort(rep(letters[1:20], times = 3)),
    Charge = rep(rep(c(2, 3), each = 15), times = 2),
    Run = rep(rep(1:3, times = 5), times = 2),
    Intensity = seq(0.1, 10, length.out = 60)
)
## With channels
to_filter_channel = data.table::data.table(
    ProteinName = "A",
    PeptideSequence = sort(rep(rep(letters[1:10], times = 3), each = 5)),
    Charge = 3,
    Run = rep(rep(1:3, each = 5), 10),
    Channel = rep(1:5, times = 30),
    Fraction = 1,
    Intensity = seq(0.1, 10, length.out = 150)
)
to_filter_channel$Intensity[1:6] = NA
# Remove feature with few measurements ----
## Nothing is missing
tinytest::expect_equal(
    MSstatsConvert:::.filterFewMeasurements(to_filter, 0, "remove", c("PeptideSequence", "Charge")),
    to_filter
)
## With a threshold:
tinytest::expect_equal(
    MSstatsConvert:::.filterFewMeasurements(to_filter, 1, "remove", c("PeptideSequence", "Charge")),
    to_filter[!(PeptideSequence %in% c("a", "b"))]
)
## With missing values
to_filter$Intensity = ifelse(to_filter$PeptideSequence %in% c("a", "b", "c") &
                                 to_filter$Run %in% 1:2, NA, to_filter$Intensity)
to_filter$Intensity[3] = NA
## Remove only features that have all missing values
tinytest::expect_equal(
    MSstatsConvert:::.filterFewMeasurements(to_filter, 0, "keep", c("PeptideSequence", "Charge")),
    to_filter[PeptideSequence != "a", ]
)
## Remove all features with less than three non-missing values
tinytest::expect_equal(
    MSstatsConvert:::.filterFewMeasurements(to_filter, 0, "remove", c("PeptideSequence", "Charge")),
    to_filter[!(PeptideSequence %in% c("a", "b", "c"))]
)
# Aggregating / summarizing features ----
to_aggregate = data.table::data.table(PeptideSequence = "A", 
                                      Run = rep(1:3, each = 2),
                                      Intensity = 1:6)
to_aggregate$Intensity[5] = NA
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultipleMeasurements(to_aggregate, max, "PeptideSequence"),
    to_aggregate[c(2, 4, 6), ]
)
to_aggregate_channel = data.table::copy(to_filter_channel)
to_aggregate_channel$Run = rep(1:6, each = 25)
to_aggregate_channel$PSM = paste(to_aggregate_channel$PeptideSequence,
                                 to_aggregate_channel$Charge,
                                 rep(1:30, each = 5), sep = "_")
to_aggregate_channel$Intensity[1:15] = 1
to_aggregate_channel$Feature = MSstatsConvert:::.combine(to_aggregate_channel$PeptideSequence,
                                                         to_aggregate_channel$Charge)
## PSM summarization works as expected
### Case 1: all PSM actually have the same intensity.
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[1:15, ], max),
    "a_3_1"
)
### Case 2: there is a score
to_aggregate_channel$Score = 1
to_aggregate_channel$Score[16:30] = rep(1:3, each = 5)
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[16:30, ], max),
    "b_3_6"
)
### Case 3: there is isolation interference
to_aggregate_channel$IsolationInterference = 1
to_aggregate_channel$IsolationInterference[31:45] = rep(1:3, each = 5)
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[31:45, ], max),
    "c_3_7"
)
### Case 4: there is IonsScore
to_aggregate_channel$IonsScore = 1
to_aggregate_channel$IonsScore[46:60] = rep(1:3, each = 5)
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[46:60, ], max),
    "d_3_12"
)
### Case 5: number of nonmissing values
to_aggregate_channel$Intensity[c(61:63, 66:67, 71)] = NA 
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[61:75, ], max),
    "e_3_15"
)
### Case 6: maximal intensity
to_aggregate_channel$Intensity[76:90] = 1:15
tinytest::expect_equal(
    MSstatsConvert:::.summarizeMultiplePSMs(to_aggregate_channel[76:90, ], max),
    "f_3_18"
)
## Aggregation works as expected
to_aggregate_channel_2 = data.table::copy(to_aggregate_channel)
to_aggregate_channel_2[, .(n_psms = data.table::uniqueN(PSM)), by = c("PeptideSequence", "Charge")]

aggregated = MSstatsConvert:::.aggregatePSMstoPeptideIons(to_aggregate_channel_2, c("PeptideSequence", "Charge"), summary_function = sum)
tinytest::expect_equal(nrow(aggregated), 5 * (data.table::uniqueN(to_aggregate_channel_2[, .(PeptideSequence, Charge, Run)])))
# Removing proteins with just a single feature ----
to_filter_channel2 = to_filter_channel
dataset_tmt2 = rbind(to_filter_channel2, 
                     data.table::data.table(ProteinName = "B", PeptideSequence = "z",
                                            Charge = 3, Run = rep(1:3, each = 5),
                                            Fraction = 1, Channel = rep(1:5, times = 3), Intensity = 1))
## Do nothing
tinytest::expect_equal(
    MSstatsConvert:::.handleSingleFeaturePerProtein(to_filter_channel2, FALSE),
    to_filter_channel2
)
## Filter when needed
tinytest::expect_equal(
    MSstatsConvert:::.handleSingleFeaturePerProtein(to_filter_channel2, TRUE),
    data.table::setkey(to_filter_channel2[!(ProteinName == "B"), ], NULL)
)

# # Utility function ----
tinytest::expect_equal(MSstatsConvert:::.combine(c("A", "B"), c("A", "B")),
                       c("A_A", "B_B"))
# General functions for cleaning features ----
## Most general function with MSstats input
tinytest::expect_equal(
    MSstatsConvert:::.cleanByFeature(to_filter, c("PeptideSequence", "Charge"),
                                     list(handle_features_with_few_measurements = "remove",
                                          summarize_multiple_psms = max)),
    to_filter[10:60, ]
)
## Most general function with MSstatsTMT function.
tinytest::expect_equal(
    MSstatsConvert:::.cleanByFeature(to_filter_channel, c("PeptideSequence", "Charge"),
                                     list(handle_features_with_few_measurements = "remove",
                                          summarize_multiple_psms = max,
                                          remove_psms_with_any_missing = TRUE)),
    to_filter_channel[PeptideSequence != "a", 
                      list(ProteinName, PeptideSequence, Charge, Run, Channel, Fraction,
                           Intensity, PSM = paste(PeptideSequence, Charge, sep = "_"))]
)
## MSstats function
tinytest::expect_equal(
    MSstatsConvert:::.cleanByFeatureMSstats(to_filter, c("PeptideSequence", "Charge"),
                                            max, "remove"),
    to_filter[10:60, ]
)
## MSstatsTMT function
tinytest::expect_equal(
    MSstatsConvert:::.cleanByFeatureTMT(to_filter_channel, c("PeptideSequence", "Charge"), max, "remove", TRUE),
    to_filter_channel[PeptideSequence != "a", 
                      list(ProteinName, PeptideSequence, Charge, Run, Channel, Fraction,
                           Intensity, PSM = paste(PeptideSequence, Charge, sep = "_"))]
)

