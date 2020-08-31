.makeBalancedDesign = function(input, fill_missing) {
    if (is.element("Channel", colnames(input))) {
        input$IsotopeLabelType = "L"
    }
    if (fill_missing) {
        all_possibilities = .getFullDesign(input)
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, list(ProteinName, feature, PeptideSequence, PrecursorCharge,
                                FragmentIon, ProductCharge, Fraction)]), 
            all.x = TRUE, by = c("feature", "Fraction"))
        all_possibilities = merge(
            all_possibilities, 
            unique(input[, list(Run, Condition, BioReplicate)]),
            all.x = TRUE, by = "Run")
        input = merge(all_possibilities, unique(input[, list(feature, Run, IsotopeLabelType, Fraction, Intensity)]),
                      all.x = TRUE, by = c("feature", "Run", "IsotopeLabelType", "Fraction"))
        # TODO: log, whether any changes were made here
    } else {
        any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, feature]))
        msg = paste("The following features have missing values in at least one run.",
                    paste(any_missing, sep = ",\n ", collapse = ",\n "))
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    if (is.element("Channel", colnames(input))) {
        input[, colnames(input) != "IsotopeLabelType", with = FALSE]
    } else {
        input
    }
}

.getFullDesign = function(input) {
    fractions = unique(input$Fraction)
    by_fraction = vector("list", length(fractions))
    for (fraction_id in seq_along(fractions)) {
        fraction = fractions[fraction_id]
        by_fraction[[fraction_id]] = data.table::as.data.table(
            expand.grid(IsotopeLabelType = unique(input[Fraction == fraction, IsotopeLabelType]),
                        feature = unique(input[Fraction == fraction, feature]),
                        Run = unique(input[Fraction == fraction, Run])))
        by_fraction[[fraction_id]]$Fraction = fraction
    }
    rbindlist(by_fraction)
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

