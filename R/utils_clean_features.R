#' Remove features with a small number of (non-missing) measurements across runs
#' @param input data.table pre-processed by one of the .cleanRaw* functions.
#' @param min_intensity minimum intensity that will be considered non-missing.
#' @param handle_few chr, if "remove", features that have less than three 
#' measurements will be removed. If "keep", only features with all missing runs
#' will be removed.
#' @param features_columns chr, vector of names of columns that define features. 
#' @return data.table
#' @keywords internal
.filterFewMeasurements = function(input, min_intensity, handle_few,
                                  feature_columns) {
    counts = input[, .(n_obs = sum(Intensity > min_intensity, 
                                   na.rm = TRUE)),
                   by = feature_columns]
    if (handle_few == "remove") {
        counts = counts[counts$n_obs > 2, ]
        getOption("MSstatsLog")("INFO", "Features with all missing measurements across runs are removed")
        getOption("MSstatsMsg")("INFO", "Features with all missing measurements across runs are removed")
    } else {
        counts = counts[counts$n_obs > 0, ]
        getOption("MSstatsLog")("INFO", "Features with 1 or two measurements across runs are removed")
        getOption("MSstatsMsg")("INFO", "Features with 1 or two measurements across runs are removed")
    }
    merge(input, counts[, feature_columns, with = FALSE],
          by = feature_columns)
}


#' Summarize multiple measurements per feature in a single run
#' @param input data.table pre-processed by one of the .cleanRaw* functions.
#' @param aggregator function that will be used to aggregate duplicated values.
#' @param feature_columns chr, vector of names of columns that define features. 
#' @return data.table
#' @keywords internal
.summarizeMultipleMeasurements = function(input, aggregator, feature_columns) {
    counts = input[, ("n_obs" = length("Intensity")), 
                   by = feature_columns]
    if(any(counts[["n_obs"]] > length(unique(input[["Run"]])))) {
        input = merge(input[, .(Intensity = aggregator(Intensity)), 
                            by = c("Run", feature_columns), with = FALSE],
                      input[, -which(colnames(input) == "Intensity"), 
                            with = FALSE],
                      by = c("Run", feature_columns)
        )
        getOption("MSstatsLog")("INFO", "Multiple measurements per run are aggregated")
        getOption("MSstatsMsg")("INFO", "Multiple measurements per run are aggregated")
    }
    input
}

#' A set of common operations for converters: remove few features and aggregate
#' 
#' This function aggregates duplicated measurements per run for all features
#' and removes features that have only missing or less than three measurements 
#' across runs. 
#' 
#' @param input
#' 
#' @param summarize_function function that will be used to aggregate multiple 
#' measurement per feature in a single run.
#' @param handle_few_measurements lgl, if TRUE, features with less than three
#' measurements across runs will be removed.
#' @return data.table
#' @keywords internal
.cleanByFeature = function(input, feature_columns, summarize_function,
                           handle_few_measurements) {   
    input = .filterFewMeasurements(input, 1, "keep", feature_columns)
    input = .summarizeMultipleMeasurements(input, summarize_function,
                                           feature_columns)
    input = .filterFewMeasurements(input, 0, handle_few_measurements,
                                   feature_columns)
    input
}


#' Remove proteins only identified by a single feature
#' @param input data.table pre-processed by one of the .cleanRaw* functions.
#' @param remove_single_feature lgl, if TRUE, proteins with a single feature
#' will be removed.
#' @param feature_columns chr, vector of names of columns that define features. 
#' @return data.table
#' @keywords internal
.handleSingleFeaturePerProtein = function(input, remove_single_feature,
                                          feature_columns) {
    counts = unique(input[, c(feature_columns, "ProteinName"), with = FALSE])
    counts = counts[, .(n_obs = .N), by = "ProteinName"]
    counts = counts[counts$n_obs > 1L, .(ProteinName)]
    if (remove_single_feature & nrow(counts) > 0) {
        input = merge(input, counts, by = "ProteinName")
        getOption("MSstatsLog")("INFO", "Proteins with a single feature are removed")
        getOption("MSstatsMsg")("INFO", "Proteins with a single feature are removed")
    }
    input
}
