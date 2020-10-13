#' Fill missing rows to create balanced design
#' @param input output of `MSstatsPreprocess`
#' @param fill_missing if TRUE, missing Intensities values will be added to data 
#' and marked as NA
#' @return data.table
#' @keywords internal
.makeBalancedDesign = function(input, fill_missing) {
    feature = NULL
    
    is_tmt = is.element("Channel", colnames(input))
    
    if (fill_missing) {
        cols = intersect(colnames(input), 
                         c("ProteinName", "feature", "PeptideSequence", "PSM", 
                           "PrecursorCharge", "FragmentIon", "ProductCharge"))
        annotation_cols = intersect(colnames(input), 
                                    c("Run", "Condition", "BioReplicate", "Channel",
                                      "Mixture", "TechRepMixture", "TechReplicate"))
        intensity_ids = intersect(c("feature", "Run", "Channel", "IsotopeLabelType",
                                    "Fraction"), colnames(input))
        if (is_tmt) {
            group_col = "Run"
            measurement_col = "Channel"
            cols = intersect(c(cols, "Fraction"),
                             colnames(input))
        } else {
            group_col = "Fraction"
            measurement_col = "Run"
        }
        all_possibilities = .getFullDesign(input, group_col, "feature",
                                           measurement_col, is_tmt)
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, c(cols, group_col), with = FALSE]), 
            all.x = TRUE, by = c("feature", group_col), allow.cartesian = TRUE)
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, annotation_cols, with = FALSE]),
            all.x = TRUE, by = unique(c("Run", measurement_col)))
        intensities = intersect(c(intensity_ids, "Intensity", "isZero"), 
                                colnames(input))
        input = merge(all_possibilities, 
                      unique(input[, intensities, with = FALSE]),
                      all.x = TRUE, by = intensity_ids)
    } else {
        if (!is_tmt) {
            any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, feature]))
            msg = paste("The following features have missing values in at least one run.",
                        paste(any_missing, sep = ",\n ", collapse = ",\n "))
            getOption("MSstatsLog")("WARN", msg)
            getOption("MSstatsMsg")("WARN", msg)
        }
    }
    input
}


#' Create a data.frame of each combination of values for given variables
#' @param input output of `MSstatsPreprocess`
#' @param group_col name of column in `input`. Combination of values of 
#' `feature_col` and `measurement_col` will be created within each unique value
#' of this column
#' @param `feature_column` name of the column that labels features
#' @param `measurement_col` name of a column with measurement labels - Runs in
#' label-free case, Channels in TMT case.
#' @param is_tmt if TRUE, data will be treated as coming from TMT experiment.
#' @importFrom data.table rbindlist
#' @return data.table
#' @keywords internal
.getFullDesign = function(input, group_col, feature_col, 
                          measurement_col, is_tmt
) {
    if (is_tmt) {
        labels = "L"
        groups = unique(input[[group_col]])
        by_group = vector("list", length(groups))
        measurements = unique(input[[measurement_col]])
        for (group_id in seq_along(groups)) {
            group = groups[group_id]
            group_filter = input[[group_col]] == group
            by_group[[group_id]] = data.table::as.data.table(
                expand.grid(labels = labels,
                            features = unique(input[[feature_col]][group_filter]),
                            measurements = measurements))
            by_group[[group_id]]$group = group
        }
        result = data.table::rbindlist(by_group)
        colnames(result) = c("IsotopeLabelType", feature_col, 
                             measurement_col, group_col)
        result[, 2:4, with = FALSE]
    } else {
        labels = unique(input[["IsotopeLabelType"]])
        groups = unique(input[[group_col]])
        by_group = vector("list", length(groups))
        measurements = unique(input[[measurement_col]])
        for (group_id in seq_along(groups)) {
            group = groups[group_id]
            group_filter = input[[group_col]] == group
            by_group[[group_id]] = data.table::as.data.table(
                expand.grid(
                    labels = labels,
                    features = unique(input[[feature_col]][group_filter]),
                    measurements = unique(input[[measurement_col]][group_filter])
                ))
            by_group[[group_id]]$group = group
        }
        result = data.table::rbindlist(by_group)
        colnames(result) = c("IsotopeLabelType", feature_col, 
                             measurement_col, group_col)
        result
    }
}


#' Get names of missing runs
#' @param input output of `MSstatsPreprocess`
#' @importFrom data.table uniqueN
#' @return data.table
#' @keywords internal
.getMissingRunsPerFeature = function(input) {
    n_measurements = Run = NULL
    
    grouping_cols = intersect(c("Fraction", "IsotopeLabelType", "feature"),
                              colnames(input))
    
    n_runs = data.table::uniqueN(input$Run)
    any_missing = input[, list(n_measurements = data.table::uniqueN(Run)),
                        by = grouping_cols]
    any_missing = any_missing[n_measurements < n_runs]
    any_missing
}

#' Check if there are duplicated measurements within run
#' @param input output of `MSstatsPreprocess`
#' @return character vector of feature labels
#' @keywords internal
.checkDuplicatedMeasurements = function(input) {
    n_measurements = feature = NULL
    counts = input[, list(n_measurements = .N),
                   by = c("Fraction", "ProteinName", "IsotopeLabelType", 
                          "feature", "Run")]
    counts = unique(counts[n_measurements > 1, as.character(feature)])
    if (length(counts) > 0) {
        msg = paste("** The following features have duplicated measurements",
                    "in some runs: please remove the duplicates.", "\n",
                    paste(counts, sep = ",\n ", collapse = ",\n "))
        # TODO: report separately for L and H
        getOption("MSstatsLog")("WARN", msg)
        stop(msg)
    }
}


#' Change labels for missing values
#' @param input output of `MSstatsPreprocess`
#' @param fix_missing missing values can be labeled by `NA`, `0` or both.
#' If `NULL`, data were processed by Skyline, so missing values will be denoted
#' by both `NA` and `0`. If "na_to_zero", `NA` values will be replaced by `0`.
#' If "zero_to_na", `0` values will be replaced by `NA`
#' @return data.table
#' @keywords internal
.fixMissingValues = function(input, fix_missing = NULL) {
    Intensity = isZero = NULL
    
    if (is.element("isZero", colnames(input))) {
        input[, Intensity := ifelse(Intensity == 0 & !is.na(Intensity), 
                                    ifelse(isZero, 0, NA), Intensity)]
        input[, Intensity := ifelse(!is.na(Intensity) & Intensity > 0 & 
                                        Intensity < 1,
                                    0, Intensity)]
    } else {
        if (!is.null(fix_missing)) {
            if (fix_missing == "zero_to_na") {
                input[, Intensity := ifelse(Intensity == 0 & !is.na(Intensity), 
                                            NA, Intensity)]
            } else {
                input[, Intensity := ifelse(is.na(Intensity), 0, Intensity)]
            }
        }
    }
    input
}
