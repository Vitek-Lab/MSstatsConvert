#' @keywords internal
.makeBalancedDesign = function(input, fill_missing) {
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
            all.x = TRUE, by = c("feature", group_col))
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, annotation_cols, with = FALSE]),
            all.x = TRUE, by = unique(c("Run", measurement_col)))
        input = merge(all_possibilities, 
                      unique(input[, c(intensity_ids, "Intensity"), with = FALSE]),
                      all.x = TRUE, by = intensity_ids)        # TODO: log, whether any changes were made here
    } else {
        any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, feature]))
        msg = paste("The following features have missing values in at least one run.",
                    paste(any_missing, sep = ",\n ", collapse = ",\n "))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    input
}


#' @keywords internal
#' @importFrom data.table rbindlist
.getFullDesign = function(input, group_col, feature_col, measurement_col, is_tmt) {
    if (is_tmt) {
        labels = "L"
    } else {
        labels = unique(input[["IsotopeLabelType"]])
    }
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
    colnames(result) = c("IsotopeLabelType", feature_col, measurement_col, group_col)
    if (is_tmt) {
        result[, 2:4, with = FALSE]
    } else {
        result
    }
}


#' @keywords internal
#' @importFrom data.table uniqueN
.getMissingRunsPerFeature = function(input) {
    n_measurements = NULL
    
    n_runs = data.table::uniqueN(input$Run)
    any_missing = input[, list(n_measurements = data.table::uniqueN(Run)),
                        by = c("Fraction", "IsotopeLabelType", "feature")]
    any_missing = any_missing[n_measurements < n_runs]
    any_missing
}


#' @keywords internal
.checkDuplicatedMeasurements = function(input) {
    n_measurements = feature = NULL
    counts = input[, list(n_measurements = .N),
                   by = c("Fraction", "ProteinName", "IsotopeLabelType", "feature", "Run")]
    counts = unique(counts[n_measurements > 1, as.character(feature)])
    if (length(counts) > 0) {
        msg = paste("The following features have duplicated measurements in some runs:",
                    "please remove the duplicates.", "\n",
                    paste(counts, sep = ",\n ", collapse = ",\n "))
        # TODO: report separately for L and H
        getOption("MSstatsLog")("WARN", msg)
        stop(msg)
    }
}
