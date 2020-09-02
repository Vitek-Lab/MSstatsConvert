#' Perform by-feature operations.
#' @param input `data.table` preprocessed by one of the cleanRaw* functions.
#' @param feature_columns character vector of names of columns that define features.
#' @param cleaning_control named list of two or three elements. See the documentation
#' for `MSstatsImport` for details.
#' @return `data.table`
#' @keywords internal 
.cleanByFeature = function(input, feature_columns, cleaning_control) {
    if (is.element("Channel", colnames(input))) {
        input = .filterFewMeasurements(
            input, 0, 
            cleaning_control[["handle_features_with_few_measurements"]],
            unique(c("PSM", "Run")))
        input = .aggregatePSMstoPeptideIons(
            input, feature_columns, 
            cleaning_control[["summarize_multiple_psms"]])
        input
    } else {
        input = .summarizeMultipleMeasurements(
            input, cleaning_control[["summarize_multiple_psms"]],
            c(feature_columns, "Run"))
        input = .filterFewMeasurements(
            input, 0, 
            cleaning_control[["handle_features_with_few_measurements"]],
            feature_columns)
    }
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
    if (is.element("Channel", colnames(input))) {
        feature_columns = c("PSM", "Run")
    }
    Intensity = n_obs = NULL
    input[, n_obs := sum(Intensity > min_intensity, na.rm = TRUE),
          by = feature_columns]
    if (handle_few == "remove") {
        cutoff = 2
        msg = "Features with one or two measurements across runs are removed"
    } else {
        cutoff = 0
        msg = "Features with all missing measurements across runs are removed"
    }
    input = input[n_obs > cutoff]
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input[, colnames(input) != "n_obs"]
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
        input[, feature := do.call(".combine", .SD), .SDcols = feature_columns]
        input[, feature_count := uniqueN(feature), by = "ProteinName"]
        input = input[feature_count > 1]
        input = input[, !(colnames(input) %in% c("feature_count", "feature")), with = FALSE]
        getOption("MSstatsLog")("INFO", "Proteins with a single feature are removed")
        getOption("MSstatsMsg")("INFO", "Proteins with a single feature are removed")
    }
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

    cols = intersect(colnames(input),
                     c("PSM", "Channel", "Intensity", "Run", "Score",
                       "IsolationInterference", "IonsScore", "n_psms"))
    input[, keep := .summarizeMultiplePSMs(.SD, summary_function), 
          by = feature_columns, .SDcols = cols]
    input = input[(PSM == keep) | is.na(keep), 
                  !(colnames(input) %in% c("keep", "feature")), 
                  with = FALSE]
    input[, n_psms := data.table::uniqueN(PSM), by = feature_columns]
    if (any(input$n_psms)) {
        input = input[, .(Intensity = mean(Intensity, na.rm = TRUE)),
                      by = setdiff(colnames(input), c("PSM", "Intensity", "n_psms",
                                                      "IsolationInterference", "Score",
                                                      "IonsScore"))]
    }
    input[, PSM := do.call(".combine", .SD), .SDcols = feature_columns]
    msg = "PSMs have been aggregated to peptide ions."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input[, colnames(input) != "n_psms", with = FALSE]
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
        is_max = sum(by_max$Intensity == max(by_max$Intensity))
        if (sum(is_max, na.rm = TRUE) == 1) {
            return(by_max$PSM[which.max(by_max$Intensity)])
        } else {
            return(NA)
        }
    }
}
