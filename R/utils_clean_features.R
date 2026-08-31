#' Perform by-feature operations.
#' @param input `data.table` preprocessed by one of the cleanRaw* functions.
#' @param feature_columns character vector of names of columns that define features.
#' @param cleaning_control named list of two or three elements. 
#' See the documentation for `MSstatsImport` for details.
#' @param anomaly_metrics character vector of quality metric column names to be 
#' used as features in an anomaly detection model.
#' @return `data.table`
#' @keywords internal 
.cleanByFeature = function(input, feature_columns, 
                           cleaning_control, anomaly_metrics = c()) {
    if (is.element("Channel", colnames(input))) {
        input = .filterFewMeasurements(
            input, 0, 
            cleaning_control[["remove_features_with_few_measurements"]],
            unique(c("PSM", "Run")))
        input = .aggregatePSMstoPeptideIons(
            input, feature_columns, 
            cleaning_control[["summarize_multiple_psms"]])
        input
    } else {
        summary_str = deparse(cleaning_control[["summarize_multiple_psms"]])
        summary_str = ifelse(grepl("sum", summary_str), "sum", "max")
        
        input = .summarizeMultipleMeasurements(
            input, cleaning_control[["summarize_multiple_psms"]],
            c(feature_columns, "Run"),
            anomaly_metrics)
        msg = paste("** Multiple measurements in a feature and a run",
                    "are summarized by summaryforMultipleRows:", summary_str)
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        input = .filterFewMeasurements(
            input, 0, 
            cleaning_control[["remove_features_with_few_measurements"]],
            feature_columns)
    }
    input
}


#' Remove features with a small number of (non-missing) measurements across runs
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param min_intensity minimum intensity that will be considered non-missing.
#' @param remove_few logical, if TRUE, features that have less than three 
#' measurements will be removed. If FALSE, only features with all missing runs
#' will be removed.
#' @param feature_columns chr, vector of names of columns that define features. 
#' @return data.table
#' @keywords internal
.filterFewMeasurements = function(input, min_intensity, remove_few,
                                  feature_columns = NULL) {
    Intensity = n_obs = NULL
    
    if (is.null(feature_columns)) {
        if (is.element("Channel", colnames(input))) {
            feature_columns = c("PSM", "Run")
        } else {
            feature_columns = intersect(colnames(input),
                                        c("PeptideModifiedSequence", "Charge",
                                          "PeptideSequence", "PrecursorCharge",
                                          "FragmentIon", "ProductCharge"))
        }
    }
    feature_columns = setdiff(feature_columns, "IsotopeLabelType")
    
    input[, n_obs := sum(Intensity > min_intensity, na.rm = TRUE),
          by = feature_columns]
    
    what = ifelse(is.element("Channel", colnames(input)),
                  "channels within each run", "runs")
    if (remove_few) {
        cutoff = 2
        msg = paste("** Features with one or two measurements across", 
                    what, "are removed.")
    } else {
        cutoff = 0
        msg = paste("** Features with all missing measurements across", 
                    what, "are removed.")
    }
    input = input[n_obs > cutoff]
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input[, colnames(input) != "n_obs", with = FALSE]
}


#' Summarize multiple measurements per feature in a single run
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param aggregator function that will be used to aggregate duplicated values.
#' @param feature_columns chr, vector of names of columns that define features.
#' @param anomaly_metrics character vector of quality metric column names 
#' to be used as features in an anomaly detection model.
#' @return `data.table`
#' @keywords internal
.summarizeMultipleMeasurements = function(input, aggregator, 
                                          feature_columns, anomaly_metrics = c()) {
    Intensity = isZero = NULL
    
    info = unique(input[, intersect(colnames(input), 
                                    c("StandardType", "ProteinName", 
                                      "PeptideModifiedSequence", "Charge",
                                      "PeptideSequence", "PrecursorCharge",
                                      "IsotopeLabelType")), 
                        with = FALSE])
    
    if (is.element("isZero", colnames(input)) & length(anomaly_metrics) == 0) {
        input = input[, list(Intensity = aggregator(Intensity, na.rm = TRUE),
                             isZero = all(isZero | is.na(Intensity)) &
                                 !all(is.na(Intensity))), 
                      by = feature_columns]
    } else if (length(anomaly_metrics) > 0){
        
        input[, row_id := .I]
        max_per_group = input[, .(max_intensity = max(Intensity, na.rm = TRUE)),
                              by = feature_columns]
        
        input = input[max_per_group, 
                         on = c(feature_columns, "Intensity" = "max_intensity"),
                         nomatch = 0L, mult = "first"]
        input[, row_id := NULL]
        
    } else {
        input = input[, list(Intensity = aggregator(Intensity, na.rm = TRUE)), 
                      by = c(feature_columns)]
    }
    merge(input, info, 
          by = intersect(colnames(input), colnames(info)), sort = FALSE)
}


#' Remove proteins only identified by a single feature
#' @param input `data.table` pre-processed by one of the .cleanRaw* functions.
#' @param remove_single_feature lgl, if TRUE, proteins with a single feature
#' will be removed.
#' @return `data.table`
#' @keywords internal
.handleSingleFeaturePerProtein = function(input, remove_single_feature) {
    feature_count = feature = NULL
    
    if (remove_single_feature) {
        feature_columns = intersect(c("PeptideSequence", "PrecursorCharge",
                                      "FragmentIon", "ProductCharge", "Charge"),
                                    colnames(input))
        input[, feature := do.call(".combine", .SD), .SDcols = feature_columns]
        input[, feature_count := uniqueN(feature), by = "ProteinName"]
        input = input[feature_count > 1]
        input = input[, !(colnames(input) %in% c("feature_count", "feature")), 
                      with = FALSE]
        getOption("MSstatsLog")("INFO", 
                                "Proteins with a single feature are removed.")
        getOption("MSstatsMsg")("INFO", 
                                "Proteins with a single feature are removed.")
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
.aggregatePSMstoPeptideIons = function(input, feature_columns, 
                                       summary_function = sum
) {
    keep = n_psms = PSM = Intensity = NULL
    
    feature_columns = unique(c(feature_columns, "Run"))
    input[, n_psms := data.table::uniqueN(PSM), by = feature_columns]
    
    if (any(input$n_psms > 1)) {
        input_duplicates = input[n_psms != 1]
        input = input[n_psms == 1]
        if (nrow(input_duplicates) > 0) {
            cols = intersect(colnames(input_duplicates),
                             c("PSM", "Channel", "Intensity", "Run", "Score",
                               "IsolationInterference", "IonsScore", "n_psms",
                               "Purity", "PeptideProphetProbability"))
            input_duplicates[, keep := .summarizeMultiplePSMs(.SD, 
                                                              summary_function),
                             by = feature_columns, .SDcols = cols]
        }
        input = rbind(input, input_duplicates, fill = TRUE)
        cols = intersect(colnames(input),
                         c("PSM", "Channel", "Intensity", "Run", "Score",
                           "IsolationInterference", "IonsScore", "n_psms",
                           "Purity", "PeptideProphetProbability"))
        input[, keep := .summarizeMultiplePSMs(.SD, summary_function), 
              by = feature_columns, .SDcols = cols]
        input = input[(PSM == keep) | is.na(keep), 
                      !(colnames(input) %in% c("keep", "feature")), 
                      with = FALSE]
        input[, n_psms := data.table::uniqueN(PSM), by = feature_columns]
        if (any(input$n_psms > 1)) {
            input = input[, list(Intensity = mean(Intensity, na.rm = TRUE)),
                          by = setdiff(colnames(input), 
                                       c("PSM", "Intensity", "n_psms",
                                         "IsolationInterference", "Score",
                                         "IonsScore", "Purity",
                                         "Purity", "PeptideProphetProbability"))]
        }
    }
    input[, PSM := do.call(".combine", .SD), 
          .SDcols = c("PeptideSequence", "PrecursorCharge")]
    msg = "** PSMs have been aggregated to peptide ions."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    input[, colnames(input) != "n_psms", with = FALSE]
}

#' Pick one PSM from a data.table of several PSMs.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param summary_function function that will be used to aggregate intensities
#' if needed.
#' @return character - label of a chosen PSM
#' @keywords internal
.summarizeMultiplePSMs = function(input, summary_function) {
    Intensity = Score = IsolationInterference = IonsScore = PSM = NULL
    
    if (all(unique(input$n_psms) == 1)) {
        return(unique(input$PSM))
    } else {
        nonmissing_counts = input[, list(n_nonmissing = sum(!is.na(Intensity))),
                                  by = c("PSM")]
        is_max = nonmissing_counts$n_nonmissing == max(nonmissing_counts$n_nonmissing, 
                                                       na.rm = TRUE)
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
        
        if ("IsolationInterference" %in% colnames(input) &
            !any(is.na(input$IsolationInterference))) {
            by_score = input[, list(score = unique(IsolationInterference)),
                             by = c("PSM")]
            is_min = by_score$score == min(by_score$score, na.rm = TRUE)
            if (sum(is_min, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.min(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_min])]
            }
        }
        
        if ("IonsScore" %in% colnames(input) & !any(is.na(input$IonsScore))) {
            by_score = input[, list(score = unique(IonsScore)),
                             by = c("PSM")]
            is_max = by_score$score == max(by_score$score, na.rm = TRUE)
            if (sum(is_max, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.max(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_max])]
            }
        }
        
        if ("Purity" %in% colnames(input) & !any(is.na(input$Purity))) {
            by_score = input[, list(score = unique(Purity)),
                             by = c("PSM")]
            is_max = by_score$score == max(by_score$score, na.rm = TRUE)
            if (sum(is_max, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.max(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_max])]
            }
        }
        
        if ("PeptideProphet.Probability" %in% colnames(input) &
            !any(is.na(input$PeptideProphet.Probability))) {
            by_score = input[, list(score = unique(PeptideProphet.Probability)),
                             by = c("PSM")]
            is_max = by_score$score == max(by_score$score, na.rm = TRUE)
            if (sum(is_max, na.rm = TRUE) == 1) {
                return(by_score$PSM[which.max(by_score$score)])
            } else {
                input = input[PSM %in% unique(by_score$PSM[is_max])]
            }
        }
        
        
        by_max = input[, list(Intensity = summary_function(Intensity, 
                                                           na.rm = TRUE)),
                       by = c("PSM")]
        is_max = by_max$Intensity == max(by_max$Intensity, na.rm = TRUE)
        if (sum(is_max, na.rm = TRUE) == 1) {
            return(by_max$PSM[which.max(by_max$Intensity)])
        } else {
            return(NA)
        }
    }
}


#' Count regex matches per element, scoring no-match and \code{NA} as 0.
#' @param x Character vector to search.
#' @param pattern Perl-compatible regex.
#' @return Integer vector of match counts, the same length as \code{x}.
#' @keywords internal
#' @noRd
.countRegexMatches = function(x, pattern) {
    lengths(regmatches(x, gregexpr(pattern, x, perl = TRUE)))
}


#' Drop peptides carrying more than one labelable residue.
#'
#' Such peptides can be partially labeled, which the two-state turnover model
#' cannot represent, so heavy and light rows are dropped together to keep the
#' light/heavy ratio unbiased.  The number of peptides removed is logged, since
#' the exclusion is otherwise invisible to the user.
#'
#' @param dt \code{data.table} with a \code{PeptideSequence} column.
#' @param residue_regex Perl-compatible regex matching one labelable residue.
#' @param strip_regex Perl-compatible regex matching label and modification
#'   annotations, removed before counting so that residue letters inside an
#'   annotation are not counted.
#' @return \code{dt} with multiply labeled rows removed.
#' @keywords internal
#' @noRd
.filterMultiplyLabeledPeptides = function(dt, residue_regex, strip_regex) {
    stripped = gsub(strip_regex, "", dt[["PeptideSequence"]], perl = TRUE)
    n_labelable = .countRegexMatches(stripped, residue_regex)
    is_multiply_labeled = n_labelable >= 2L

    if (any(is_multiply_labeled)) {
        # Count distinct peptides on the stripped sequence, so that the heavy
        # and light forms of one peptide are not reported as two.
        msg = paste("**", data.table::uniqueN(stripped[is_multiply_labeled]),
                    "peptide(s) with more than one labelable residue were",
                    paste0("removed (", sum(is_multiply_labeled), " row(s))."),
                    "Turnover analysis is currently limited to peptides with",
                    "exactly one labelable residue.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }

    dt[!is_multiply_labeled, ]
}


#' Classify IsotopeLabelType from peptide sequence patterns.
#'
#' Shared by the Spectronaut and DIANN turnover workflows.  Each peptide is
#' classified as heavy (\code{"H"}), light (\code{"L"}), or unlabeled
#' (\code{NA}) by matching regexes against \code{PeptideSequence}.
#'
#' Spectronaut mode (\code{labeled_aa_regex}) strips bracket modifications
#' first, then infers light from a bare labelable residue with no heavy tag.
#' DIANN mode (\code{light_regex}) matches heavy and light patterns directly
#' against the modified sequence.
#'
#' @param dt \code{data.table} with a \code{PeptideSequence} column.
#' @param heavy_regex Perl-compatible regex matching heavy-labeled sequences.
#' @param light_regex Perl-compatible regex explicitly matching light-labeled
#'   sequences.  Required for DIANN mode; must be non-\code{NULL}.
#' @param labeled_aa_regex Perl-compatible regex for bare labeled amino acids
#'   applied to the bracket-stripped sequence to detect light peptides.
#'   Required for Spectronaut mode; must be non-\code{NULL}.
#'   Exactly one of \code{light_regex} and \code{labeled_aa_regex} must be
#'   supplied.
#' @return \code{dt} with \code{IsotopeLabelType} column added or updated.
#' @keywords internal
#' @noRd
.classifyIsotopeLabelType = function(dt, heavy_regex,
                                      light_regex = NULL,
                                      labeled_aa_regex = NULL) {
    IsotopeLabelType = PeptideSequence = StrippedSequence = NULL

    if (!is.null(labeled_aa_regex)) {
        dt[, StrippedSequence := gsub("\\[.*?\\]", "", PeptideSequence)]
        dt[, IsotopeLabelType := data.table::fcase(
            grepl(heavy_regex, PeptideSequence, perl = TRUE), "H",
            grepl(labeled_aa_regex, StrippedSequence, perl = TRUE), "L",
            default = NA_character_
        )]
        dt[, StrippedSequence := NULL]
    } else if (!is.null(light_regex)) {
        dt[, IsotopeLabelType := data.table::fcase(
            grepl(heavy_regex, PeptideSequence, perl = TRUE), "H",
            grepl(light_regex, PeptideSequence, perl = TRUE), "L",
            default = NA_character_
        )]
    }

    dt
}


#' Fix invalid intensities: infinite to NA, between 0 and 1 to 0
#' @param input data.table
#' @return data.table
#' @keywords internal
.adjustIntensities = function(input) {
    Intensity = isZero = NULL
    
    if (is.element("isZero", colnames(input))) {
        input[, isZero := ifelse(Intensity > 0 & Intensity <= 1, TRUE, isZero)]
    }
    input[, Intensity := ifelse(is.finite(Intensity), Intensity, NA)]
    input[, Intensity := ifelse(Intensity > 0 & Intensity <= 1, 0, Intensity)]
}
