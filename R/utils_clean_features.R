#' Perform by-feature operations.
#' @param input `data.table` preprocessed by one of the cleanRaw* functions.
#' @param feature_columns character vector of names of columns that define features.
#' @param cleaning_control named list of two or three elements. See the documentation
#' for `MSstatsImport` for details.
#' @return `data.table`
#' @keywords internal 
.cleanByFeature = function(input, feature_columns, cleaning_control) {
    if (is.element("Channel", colnames(input))) {
        .cleanByFeatureTMT(input, feature_columns, 
                           cleaning_control[["summarize_multiple_psms"]],
                           cleaning_control[["handle_features_with_few_measurements"]])
    } else {
        .cleanByFeatureMSstats(input, feature_columns,
                               cleaning_control[["summarize_multiple_psms"]],
                               cleaning_control[["handle_features_with_few_measurements"]])
    }
}


#' A set of common operations for converters: remove few features and aggregate
#' 
#' This function aggregates duplicated measurements per run for all features
#' and removes features that have only missing or less than three measurements 
#' across runs. 
#' 
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param summary_function function that will be used to aggregate multiple 
#' measurement per feature in a single run.
#' @param handle_few_measurements lgl, if TRUE, features with less than three
#' measurements across runs will be removed.
#' @return `data.table`
#' @keywords internal
.cleanByFeatureMSstats = function(input, feature_columns, summary_function,
                                  handle_few_measurements) {
    input = .summarizeMultipleMeasurements(input, summary_function,
                                           c(feature_columns, "Run"))
    input = .filterFewMeasurements(input, 0, handle_few_measurements,
                                   feature_columns)
    input
}


#' Remove features with a small number of (non-missing) measurements across runs
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param min_intensity minimum intensity that will be considered non-missing.
#' @param handle_few chr, if "remove", features that have less than three 
#' measurements will be removed. If "keep", only features with all missing runs
#' will be removed.
#' @param features_columns chr, vector of names of columns that define features. 
#' @return data.table
#' @keywords internal
.filterFewMeasurements = function(input, min_intensity, handle_few,
                                  feature_columns) {
    Intensity = n_obs = NULL
    counts = input[, list(n_obs = sum(Intensity > min_intensity, 
                                      na.rm = TRUE)),
                   by = feature_columns]
    if (handle_few == "remove") {
        counts = counts[n_obs > 2, ]
        msg = "Features with one or two measurements across runs are removed"
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    } else {
        counts = counts[n_obs > 0, ]
        msg = "Features with all missing measurements across runs are removed"
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    merge(input, unique(counts[, feature_columns, with = FALSE]),
          by = feature_columns, sort = FALSE)
}


#' Summarize multiple measurements per feature in a single run
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param aggregator function that will be used to aggregate duplicated values.
#' @param feature_columns chr, vector of names of columns that define features. 
#' @return `data.table`
#' @keywords internal
.summarizeMultipleMeasurements = function(input, aggregator, feature_columns) {
    Intensity = NULL
    
    info = unique(input[, intersect(colnames(input), 
                                    c("StandardType", "ProteinName", 
                                      "PeptideModifiedSequence", "Charge",
                                      "PeptideSequence", "PrecursorCharge",
                                      "IsotopeLabelType")), 
                        with = FALSE])
    input = input[, list(Intensity = aggregator(Intensity, na.rm = TRUE)), 
                  by = feature_columns]
    merge(input, info, by = intersect(colnames(input), colnames(info)), sort = FALSE)
}


#' Perform by-feature operations for TMT data.
#' @inheritParams .cleanByFeatureMSstats
#' @param remove_any_missing If TRUE, features that have any missing values in a
#' run will be removed from that run.
#' @return `data.table`
#' @keywords internal
.cleanByFeatureTMT = function(input, feature_columns, summary_function,
                              remove_few, remove_any_missing) {
    PSM = PeptideSequence = PrecursorCharge = NULL
    
    input = .filterFewMeasurements(input, 1, remove_few,
                                   unique(c("PSM", "Run")))
    input = .aggregatePSMstoPeptideIons(input, feature_columns, summary_function)
    input[, PSM := paste(PeptideSequence, PrecursorCharge, sep = "_")]
    input
}


#' Remove proteins only identified by a single feature
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param remove_single_feature lgl, if TRUE, proteins with a single feature
#' will be removed.
#' @return `data.table`
#' @keywords internal
.handleSingleFeaturePerProtein = function(input, remove_single_feature) {
    if (remove_single_feature) {
        feature_columns = intersect(c("PeptideSequence", "PrecursorCharge",
                                      "FragmentIon", "ProductCharge", "Charge"),
                                    colnames(input))
        input[, feature := do.call(".combine", feature_columns)]
        input[, feature_count := uniqueN(feature), by = "ProteinName"]
        input = input[feature_count > 1]
        input = input[, !(colnames(input) %in% c("feature_count", "feature")), with = FALSE]
        getOption("MSstatsLog")("INFO", "Proteins with a single feature are removed")
        getOption("MSstatsMsg")("INFO", "Proteins with a single feature are removed")
    }
    input
}

#' Remove features for which all channels are missing in one run.
#' @param input data.table preprocessed by one of the cleanRaw* functions.
#' @param feature_columns chr, vector of names of columns that define features.
#' @return `data.table`
#' @keywords internal
.removeMissingAllChannels = function(input, feature_columns) {
    Intensity = NotAllMissing = NULL
    cols = unique(c(feature_columns, "Run"))
    input[, NotAllMissing := !(all(is.na(Intensity)) | all(Intensity == 0)),
          by = cols]
    input = input[(NotAllMissing)]
    input = input[, colnames(input) != "NotAllMissing", with = FALSE]
    msg = "PSMs, that have all zero intensities across channels in each run, are removed."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input
}


#' @keywords internal
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
.aggregatePSMstoPeptideIons = function(input, feature_columns, summary_function = sum) {
    feature = keep = n_psms = PSM = NULL
    
    feature_columns = unique(c(feature_columns, "Run"))
    input[, n_psms := data.table::uniqueN(PSM), by = feature_columns]
    input[, feature := do.call(".combine", .SD), 
          .SDcols = feature_columns]
    
    cols = intersect(colnames(input),
                     c("feature", "PSM", "Channel", "Intensity", "Run", "Score",
                       "IsolationInterference", "IonsScore", "n_psms"))
    input[, keep := .summarizeMultiplePSMs(.SD, summary_function), 
          by = feature_columns, .SDcols = cols]
    input = input[PSM == keep, 
                  !(colnames(input) %in% c("keep", "feature")), 
                  with = FALSE]
    
    msg = "PSMs have been aggregated to peptide ions."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input
}

#' Pick one PSM from a data.table of several PSMs.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param summary_function function that will be used to aggregate intensities
#' if needed.
#' @keywords internal
.summarizeMultiplePSMs = function(input, summary_function) {
    Intensity = Score = IsolationInterference = IonsScore = NULL
    
    if (all(unique(input$n_psms) == 1)) {
        return(unique(input$PSM))
    } else {
        unique_counts = nrow(unique(input[, colnames(input) != "PSM", with = FALSE]))
        if (unique_counts == (nrow(input) / unique(input$n_psms))) {
            return(input$PSM[1])
        }
        
        nonmissing_counts = input[, list(n_nonmissing = sum(!is.na(Intensity))),
                                  by = c("PSM")]
        is_max = nonmissing_counts$n_nonmissing == max(nonmissing_counts$n_nonmissing, na.rm = TRUE)
        if (sum(is_max, na.rm = TRUE) == 1) {
            return(nonmissing_counts$PSM[which.max(nonmissing_counts$n_nonmissing)])
        } else {
            input = input[PSM %in% unique(nonmissing_counts$PSM[is_max])]
        }
        
        if ("Score" %in% colnames(input)) {
            by_score = input[, list(score = unique(Score)),
                             by = c("PSM")]
            is_max = by_score$score == max(by_score$score, na.rm = TRUE)
            if (sum(is_max, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.max(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_max])]
            }
        }
        
        if ("IsolationInterference" %in% colnames(input)) {
            by_score = input[, list(score = unique(IsolationInterference)),
                             by = c("PSM")]
            is_min = by_score$score == min(by_score$score, na.rm = TRUE)
            if (sum(is_min, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.min(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_min])]
            }
        }
        
        if ("IonsScore" %in% colnames(input)) {
            by_score = input[, list(score = unique(IonsScore)),
                             by = c("PSM")]
            is_max = sum(by_score$score == max(by_score$score, na.rm = TRUE))
            if (sum(is_max, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.max(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_max])]
            }
        }
        
        by_max = input[, list(Intensity = summary_function(Intensity, na.rm = TRUE)),
                       by = c("PSM")]
        return(by_max$PSM[which.max(by_max$Intensity)])
    }
}
