#' Check if there are overlapping features and remove if needed
#' @param input data.table preprocessed by one of the .cleanRaw* functions and
#' merged with annotation.
#' @return `data.table`
#' @keywords internal
.handleFractions = function(input) {
    fractions = unique(input$Fraction)
    if (length(fractions) > 1 & is.element("Channel", colnames(input))) {
        input = .combineFractions(input)
        msg = "Fractions belonging to same mixture have been combined."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    input
}

#' Remove peptide ions overlapped among multiple fractions of the same biological mixture
#' @inheritParams .handleFractions
#' @return `data.table`
#' @keywords internal
.combineFractions = function(input) {
    techrun = feature = id = Mixture = TechRepMixture = Run = NULL
    
    cols = intersect(colnames(input), c("ProteinName", "PeptideSequence", 
                                        "Charge", "PrecursorCharge"))
    input[, techrun := paste(Mixture, TechRepMixture, sep = "_")]
    input[, feature := do.call(".combine", .SD), 
          .SDcols = cols]
    input[, id := paste(feature, Run, sep = "_")]
    
    unoverlapped_list = vector("list", length(unique(input$techrun)))
    names(unoverlapped_list) = unique(input$techrun)
    for (technical_run in unique(input$techrun)) {
        single_run = input[techrun == technical_run & !is.na(Intensity), ]
        features_to_remove = .getOverlappingFeatures(single_run)
        if (length(features_to_remove) > 0) {
            single_run = .filterOverlapped(single_run, mean, features_to_remove)
            features_to_remove = .getOverlappingFeatures(single_run)
            msg = paste("For peptides overlapped between fractions of",
                        technical_run, "use the fraction with maximal average abundance.")
            getOption("MSstatsLog")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg)
            if (length(features_to_remove) > 0) {
                single_run = .filterOverlapped(single_run, sum, features_to_remove)
                features_to_remove = .getOverlappingFeatures(single_run)
                msg = paste("For peptides overlapped between fractions of",
                            technical_run, "use the fraction with maximal summation abundance.")
                getOption("MSstatsLog")("INFO", msg)
                getOption("MSstatsMsg")("INFO", msg)
                if (length(features_to_remove) > 0) {
                    single_run = .filterOverlapped(single_run, max, features_to_remove)
                    msg = paste("For peptides overlapped between fractions of",
                                technical_run, "use the fraction with maximal abundance.")
                    getOption("MSstatsLog")("INFO", msg)
                    getOption("MSstatsMsg")("INFO", msg)
                    if (data.table::uniqueN(input$Run) > 1) {
                        single_run = single_run[, list(Intensity = mean(Intensity, na.rm = TRUE)),
                                                by = setdiff(colnames(single_run), 
                                                             c("Run", "Intensity",
                                                               "Fraction", "id", "n_psms",
                                                               "QuanInfo", "IonsScore",
                                                               "IsolationInterference"))]
                    }
                }
            }
        }
        unoverlapped_list[[technical_run]] = single_run[, !(colnames(single_run) %in% c("feature", "run", "Fraction", "techrun", "id")),
                                                        with = FALSE]
    }
    input = rbindlist(unoverlapped_list)
    input$Run = paste(input$Mixture, input$TechRepMixture, sep = "_")
    input
}


#' Get features that are overlapped among multiple runs
#' @inheritParams .handleFractions
#' @return `data.table`
#' @keywords internal
.getOverlappingFeatures = function(input) {
    Run = feature = n_runs = NULL
    
    count_fractions = input[, list(n_runs = uniqueN(Run)),
                            by = "feature"]
    count_fractions[n_runs > 1, feature]
}


#' Remove overlapped features
#' @inheritParams .handleFractions
#' @param summary_function summary function (mean, sum, max) that will be used
#' to pick one feature from multiple overlapping features
#' @param overlapped_features features that overlap.
#' @return `data.table`
#' @keywords internal
.filterOverlapped = function(input, summary_function, overlapped_features) {
    Intensity = id = agg_intensity = max_intensity = feature =  NULL
    
    overlapped = input[feature %in% overlapped_features]
    overlapped[, agg_intensity := summary_function(Intensity, na.rm = TRUE),
               by = c("feature", "id")]
    overlapped[, max_intensity := max(agg_intensity),
               by = "feature"]
    overlapped = overlapped[agg_intensity != max_intensity]
    input[!(id %in% overlapped$id), 
          !(colnames(input) %in% c("agg_intensity", "max_intensity")), 
          with = FALSE]
}
