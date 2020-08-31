.makeBalancedDesign = function(input, fill_missing) {
    is_tmt = is.element("Channel", colnames(input))
    if (is_tmt) {
        input$IsotopeLabelType = "L"
    } else {
        input$Channel = 1
    }
    if (fill_missing) {
        cols = intersect(colnames(input), 
                         c("ProteinName", "feature", "PeptideSequence", 
                           "PrecursorCharge", "FragmentIon", "ProductCharge", 
                           "Fraction"))
        annotation_cols = intersect(colnames(input), 
                                    c("Run", "Condition", "BioReplicate", "Channel",
                                      "Mixture", "TechRepMixture", "TechReplicate"))
        all_possibilities = .getFullDesign(input)
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, cols, with = FALSE]), 
            all.x = TRUE, by = c("feature", "Fraction"))
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, annotation_cols, with = FALSE]),
            all.x = TRUE, by = c("Run", "Channel"))
        input = merge(all_possibilities, 
                      unique(input[, list(feature, Run, Channel, IsotopeLabelType, Fraction, Intensity)]),
                      all.x = TRUE, by = c("feature", "Run", "Channel", "IsotopeLabelType", "Fraction"))
        # TODO: log, whether any changes were made here
    } else {
        any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, feature]))
        msg = paste("The following features have missing values in at least one run.",
                    paste(any_missing, sep = ",\n ", collapse = ",\n "))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    if (is_tmt) {
        input[, colnames(input) != "IsotopeLabelType", with = FALSE]
    } else {
        input[, colnames(input) != "Channel", with = FALSE]
    }
}

.getFullDesign = function(input) {
    fractions = unique(input$Fraction)
    by_fraction = vector("list", length(fractions))
    channels = unique(input$Channel)
    for (fraction_id in seq_along(fractions)) {
        fraction = fractions[fraction_id]
        by_fraction[[fraction_id]] = data.table::as.data.table(
            expand.grid(IsotopeLabelType = unique(input[, IsotopeLabelType]),
                        feature = unique(input[Fraction == fraction, feature]),
                        Run = unique(input[Fraction == fraction, Run]),
                        Channel = channels))
        by_fraction[[fraction_id]]$Fraction = fraction
    }
    data.table::rbindlist(by_fraction)
}

.getMissingRunsPerFeature = function(input) {
    n_runs = data.table::uniqueN(input$Run)
    any_missing = input[, list(n_measurements = data.table::uniqueN(Run)),
                        by = c("Fraction", "IsotopeLabelType", "feature")]
    any_missing = any_missing[n_measurements < n_runs]
    any_missing
}


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

