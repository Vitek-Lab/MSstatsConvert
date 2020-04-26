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

#' Remove features for which all channels are missing in one run.
#' @param input data.table preprocessed by one of the cleanRaw* functions.
#' @param feature_columns chr, vector of names of columns that define features.
#' @return data.table
#' @keywords internal
.removeMissingAllChannels = function(input, feature_columns) {
    cols = c("ProteinName", feature_columns, "Run")
    non_missing = input[, .(AllMissing = all(is.na(Intensity)) | all(Intensity == 0)),
                        by = cols]
    non_missing = non_missing[!non_missing$AllMissing, cols, with = FALSE]
    input = merge(input, non_missing, by = cols)
    msg = "PSMs, that have all zero intensities across channels in each run, are removed."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input
}


.removeAnyMissingInRun = function(input, feature_columns, remove) {
    # TODO: probably this code will be more general with all(!is.na(Intensity))
    if (remove) {
        cols = c("ProteinName", feature_columns, "Run")
        n_channels = length(unique(input$Channel))
        counts = input[, .(n_not_missing = sum(!is.na(Intensity))),
                       by = cols]
        counts = counts[counts$n_not_missing == n_channels, ]
        input = merge(input, counts, by = cols)
        msg = "Rows which has any missing value within a run were removed from that run."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    input
}


.cleanByFeatureTMT = function(input, feature_columns, summary_function,
                              remove_few, remove_any_missing) {
    input = .removeAnyMissingInRun(input, feature_columns, remove_any_missing)
    # TODO: MSstats had a threshold on Intensity - is it not necessary here?
    input = .filterFewMeasurements(input, 0, remove_few,
                                   c("ProteinName", feature_columns, "Run"))
    input = .aggregatePSMstoPeptideIons(input, feature_columns, summary_function)
    input
}

.combine = function(...) {
    paste(..., sep = "_")  
} 

#' Aggregate multiple PSMs to a single peptide ion.
#' @param input data.table preprocessed by one of the cleanRaw* functions.
#' @param feature_columns chr, names of columns that define features.
#' @param summary_function function that will be used to aggregate intensities
#' if needed.
#' @return data.table
#' @keywords internal
.aggregatePSMstoPeptideIons <- function(input, feature_columns, summary_function = sum) {
    feature_columns = cols = c("ProteinName", feature_columns, "Run")
    n_channels = length(unique(input$Channel))
    cols = intersect(colnames(duplicated),
                     c("Feature", "PSM", "Channel", "Intensity", "Run", "Score",
                       "IsolationInterference", "IonsScore"))
    
    counts = input[, .(n_measurements = .N),
                   by = feature_columns]
    duplicated = counts[counts$n_measurements > n_channels, ]
    not_duplicated = merge(input, counts[counts$n_measurements == n_channels, 
                                         feature_columns, with = FALSE],
                           by = feature_columns)
    duplicated = merge(input, duplicated[, feature_columns, with = FALSE], 
                       by = feature_columns)
    duplicated[, Feature := do.call(".combine", .SD), 
               .SDcols = feature_columns]
    duplicated[, keep := .summarizeMultiplePSMs(.SD, summary_function), 
               by = "Feature", .SDcols = cols]
    duplicated = duplicated[duplicated$PSM == duplicated$keep, 
                            colnames(duplicated) != "keep"]
    input = rbind(duplicated, not_duplicated)
    
    msg = "PSMs have been aggregated to peptide ions."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input
}

#' Pick one PSM from a data.table of several PSMs.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param summary_function function that will be used to aggregate intensities
#' if needed.
.summarizeMultiplePSMs = function(input, summary_function) {
    unique_counts = input[, .(n_unique = uniqueN(Intensity)),
                          by = c("Feature", "Run", "Channel")]
    if (all(unique_counts$n_unique == 1L)) {
        return(input$PSM[1])
    }
    
    if ("Score" %in% colnames(input)) {
        by_score = input[, .(score = unique(Score)),
                         by = c("Feature", "Run", "PSM")]
        if (uniqueN(by_score$score) != 1) {
            return(by_score$PSM[which.max(by_score$score)])
        }
    }
    
    if ("IsolationInterference" %in% colnames(input)) {
        by_score = input[, .(score = unique(IsolationInterference)),
                         by = c("Feature", "Run", "PSM")]
        if (uniqueN(by_score$score) != 1) {
            return(by_score$PSM[which.min(by_score$score)])
        }
    }
    
    if ("IonsScore" %in% colnames(input)) {
        by_score = input[, .(score = unique(IonsScore)),
                         by = c("Feature", "Run", "PSM")]
        if (uniqueN(by_score$score) != 1) {
            return(by_score$PSM[which.max(by_score$score)])
        }
    }
    
    nonmissing_counts = input[, .(n_nonmissing = sum(!is.na(Intensity))),
                              by = c("Feature", "Run", "PSM")]
    if (uniqueN(nonmissing_counts$n_nonmissing) != 1) {
        return(nonmissing_counts$PSM[which.max(nonmissing_counts$n_nonmissing)])
    }
    
    by_max = input[, .(Intensity = summary_function(Intensity, na.rm = TRUE)),
                   by = c("Feature", "Run", "PSM")]
    return(by_max$PSM[which.max(by_max$Intensity)])
}
