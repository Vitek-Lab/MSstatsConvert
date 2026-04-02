# Balanced design
filter_NA = function(dt) {
    dt[!is.na(Intensity)]
}
## Label-free
### Example data
n_runs = 5
n_labels = 2
n_features = 5
test_data_1 = data.table::data.table(
    ProteinName = "A",
    feature = rep(letters[1:n_features], each = n_runs * n_labels),
    Run = rep(rep(1:n_runs, times = n_features), each = n_labels),
    IsotopeLabelType = rep(c(c("L", "H")[1:n_labels]), times = n_features * n_runs),
    Fraction = 1,
    Intensity = runif(n_runs * n_labels * n_features)
)
### Some rows are missing
test1_na = data.table::copy(test_data_1)[order(ProteinName, feature, 
                                               IsotopeLabelType, Run, Fraction)]
test1_na[, Intensity := ifelse(feature == "a" & Run %in% c(1, 3, 4) & 
                                   IsotopeLabelType == "L", 
                               NA, Intensity)]
test1_nona = filter_NA(test1_na)
expect_equal(
    MSstatsConvert:::.makeBalancedDesign(test1_nona, TRUE)[order(ProteinName, feature, 
                                                                 IsotopeLabelType, Run, Fraction), Intensity],
    test1_na$Intensity
)
### List of missing values is returned when fill_missing = FALSE
expect_stdout(MSstatsConvert:::.makeBalancedDesign(test1_nona, FALSE))
expect_equal(MSstatsConvert:::.getMissingRunsPerFeature(test1_nona)$feature,
             "a")
### All rows in one label are missing
test2_na = data.table::copy(test_data_1)[order(ProteinName, feature, 
                                               IsotopeLabelType, Run, Fraction)]
test2_na[, Intensity := ifelse(feature == "a" & IsotopeLabelType == "L", 
                               NA, Intensity)]
test2_nona = filter_NA(test2_na)
expect_equal(
    MSstatsConvert:::.makeBalancedDesign(test2_nona, TRUE)[order(ProteinName, feature, 
                                                                 IsotopeLabelType, Run, Fraction), Intensity],
    test2_na$Intensity
)
## Labeled
n_labels = 4
test_data_tmt = data.table::data.table(
    ProteinName = "A",
    feature = rep(letters[1:n_features], each = n_runs * n_labels),
    Run = rep(rep(1:n_runs, times = n_features), each = n_labels),
    Channel = rep(c(letters[1:n_labels]), times = n_features * n_runs),
    Intensity = runif(n_runs * n_labels * n_features)
)
test_data_tmt_na = data.table::copy(test_data_tmt)[order(ProteinName, feature, 
                                                         Run, Channel)]
set.seed(100)
test_data_tmt_na$Intensity[sample(1:nrow(test_data_tmt), 15)] = NA
test_data_tmt_nona = filter_NA(test_data_tmt_na)
expect_equal(
    MSstatsConvert:::.makeBalancedDesign(test_data_tmt_nona, TRUE)[order(ProteinName, feature, 
                                                                         Run, Channel), Intensity],
    test_data_tmt_na$Intensity
)
# Duplicated measurements
## No error when there are no duplicates
no_duplicates = data.table::data.table(
    ProteinName = "A",
    IsotopeLabelType = "L",
    Fraction = 1,
    feature = "a",
    Run = 1:6,
    Intensity = 1:6
)
expect_silent(MSstatsConvert:::.checkDuplicatedMeasurements(no_duplicates))
expect_error(MSstatsConvert:::.checkDuplicatedMeasurements(
    rbind(no_duplicates,
          no_duplicates[6, ])
))
# Missing values
no_duplicates2 = data.table::copy(no_duplicates)
no_duplicates2$Intensity[c(1, 6)] = 0
no_duplicates3 = data.table::copy(no_duplicates)
no_duplicates3$Intensity[c(1, 6)] = NA
no_duplicates4 = data.table::copy(no_duplicates)
no_duplicates4$Intensity[c(1, 6)] = 0
no_duplicates4$isZero = c(TRUE, rep(FALSE, 5))
## Zeros are converted to NA
expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates2), "zero_to_na"),
             no_duplicates3)
## NAs are converted to 0
expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates3), "na_to_zero"),
             no_duplicates2)
## For Skyline, 0 that are a result of sum(NA, na.rm = T) are replaced by NA
expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates4))$Intensity,
             c(0, 2:5, NA))
## Do nothing if no isZero column and fix_missing = NULL
expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates3), NULL)$Intensity,
             no_duplicates3$Intensity)

# .getFullDesign tests

# Test 1: H/L features expand over all runs in the group even if a label was only
# observed in a subset of runs. NA-labeled features also expand to all group runs.
# "b" has H only in run 1 and L only in run 2; "a" (NA) only in run 3.
# All three should be filled out to all runs {1, 2, 3} in the group.
gfd_mixed = data.table::data.table(
    feature          = c("b", "b", "a"),
    Run              = c(  1,   2,   3),
    IsotopeLabelType = c("H", "L",  NA),
    Fraction         = 1L
)
gfd_result = MSstatsConvert:::.getFullDesign(
    gfd_mixed, group_col = "Fraction", feature_col = "feature",
    measurement_col = "Run", is_tmt = FALSE)
a_rows = gfd_result[gfd_result$feature == "a", ]
expect_true(all(is.na(a_rows$IsotopeLabelType)))
expect_equal(nrow(a_rows), 3L)
b_rows = gfd_result[gfd_result$feature == "b", ]
expect_equal(nrow(b_rows), 6L)
expect_true(all(sort(unique(b_rows$IsotopeLabelType)) == c("H", "L")))

# Test 2: All features have NA IsotopeLabelType — full design contains only NA labels,
# not the empty H/L set that na.omit() would strip from the labels vector.
gfd_all_na = data.table::data.table(
    feature          = c("x", "x", "y", "y"),
    Run              = c(  1,   2,   1,   2),
    IsotopeLabelType = NA_character_,
    Fraction         = 1L
)
gfd_na_result = MSstatsConvert:::.getFullDesign(
    gfd_all_na, group_col = "Fraction", feature_col = "feature",
    measurement_col = "Run", is_tmt = FALSE)
expect_true(all(is.na(gfd_na_result$IsotopeLabelType)))
expect_equal(nrow(gfd_na_result), 4L)

# Test 3: NA-labeled features expand over ALL runs in the group, not just the runs
# where they were observed. Feature "p" appears in runs 1 and 2 with H/L;
# feature "q" appears only in run 3 with NA, but should be filled out to runs 1-3.
gfd_separate_runs = data.table::data.table(
    feature          = c("p", "p", "p", "p", "q"),
    Run              = c(  1,   1,   2,   2,   3),
    IsotopeLabelType = c("L", "H", "L", "H",  NA),
    Fraction         = 1L
)
gfd_sep_result = MSstatsConvert:::.getFullDesign(
    gfd_separate_runs, group_col = "Fraction", feature_col = "feature",
    measurement_col = "Run", is_tmt = FALSE)
p_rows = gfd_sep_result[gfd_sep_result$feature == "p", ]
expect_equal(nrow(p_rows), 6)
q_rows = gfd_sep_result[gfd_sep_result$feature == "q", ]
expect_equal(nrow(q_rows), 3L)
expect_true(all(is.na(q_rows$IsotopeLabelType)))
expect_equal(sort(q_rows$Run), c(1L, 2L, 3L))

