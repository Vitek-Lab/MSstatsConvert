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
tinytest::expect_equal(
    MSstatsConvert:::.makeBalancedDesign(test1_nona, TRUE)[order(ProteinName, feature, 
                                                                 IsotopeLabelType, Run, Fraction), Intensity],
    test1_na$Intensity
)
### All rows in one label are missing
test2_na = data.table::copy(test_data_1)[order(ProteinName, feature, 
                                               IsotopeLabelType, Run, Fraction)]
test2_na[, Intensity := ifelse(feature == "a" & IsotopeLabelType == "L", 
                               NA, Intensity)]
test2_nona = filter_NA(test2_na)
tinytest::expect_equal(
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
tinytest::expect_equal(
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
tinytest::expect_silent(MSstatsConvert:::.checkDuplicatedMeasurements(no_duplicates))
tinytest::expect_error(MSstatsConvert:::.checkDuplicatedMeasurements(
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
tinytest::expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates2), "zero_to_na"),
                       no_duplicates3)
## NAs are converted to 0
tinytest::expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates3), "na_to_zero"),
                       no_duplicates2)
## For Skyline, 0 that are a result of sum(NA, na.rm = T) are replaced by NA
tinytest::expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates4))$Intensity,
                       c(0, 2:5, NA))
## Do nothing if no isZero column and fix_missing = NULL
tinytest::expect_equal(MSstatsConvert:::.fixMissingValues(data.table::copy(no_duplicates3), NULL)$Intensity,
                       no_duplicates3$Intensity)
