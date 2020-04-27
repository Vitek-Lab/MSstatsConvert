#' Check if there are overlapping features and remove if needed
#' @param input data.table preprocessed by one of the .cleanRaw* functions and
#' merged with annotation.
#' @return data.table
#' @keywords internal
.handleFractions = function(input) {
    fractions = unique(input$Fraction)
    if (length(fractions) > 1) {
        input = .combineFractions(input)
        msg = "Fractions belonging to same mixture have been combined."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    input
}

#' Remove peptide ions overlapped among multiple fractions of the same biological mixture
#' @inheritParams .handleFractions
#' @return data.table
#' @keywords internal
.combineFractions = function(input) {
    input[, techrun := paste(Mixture, TechRepMixture, sep = "_")]
    input[, feature := do.call(".combine", .SD), 
          .SDcols = c("ProteinName", "PeptideSequence", "Charge")]
    input[, id := paste(feature, Run, sep = "_")]
    
    unoverlapped_list = vector("list", length(unique(input$techrun)))
    names(unoverlapped_list) = unique(input$techrun)
    for (technical_run in unique(input$techrun)) {
        single_run = input[input$techrun == technical_run & !is.na(input$Intensity), ]
        features_to_remove = .getOverlappingFeatures(single_run)
        if (nrow(features_to_remove) > 0) {
            single_run = .filterOverlapped(single_run, mean, features_to_remove)
            features_to_remove = .getOverlappingFeatures(single_run)
            if (nrow(features_to_remove) > 0) {
                single_run = .filterOverlapped(single_run, sum, features_to_remove)
                features_to_remove = .getOverlappingFeatures(single_run)
                if (nrow(features_to_remove) > 0) {
                    single_run = .filterOverlapped(single_run, max, features_to_remove)
                }
            }
        }
        unoverlapped_list[[technical_run]] = single_run[, !(colnames %in% c("feature", "techrun", "id")),
                                                        with = FALSE]
        msg = paste("For peptides overlapped between fractions of",
                    technical_run, "use the fraction with maximal average abundance.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    rbindlist(unoverlapped_list)
}


#' Get features that are overlapped among multiple runs
#' @inheritParams .handleFractions
#' @return data.table
#' @keywords internal
.getOverlappingFeatures = function(input) {
    count_fractions = input[, .(n_runs = uniqueN(Run)),
                            by = "feature"]
    count_fractions[count_fractions$n_runs > 1, .(feature)]
}


#' Remove overlapped features
#' @inheritParams .handleFraction
#' @param summary_function summary function (mean, sum, max) that will be used
#' to pick one feature from multiple overlapping features
#' @param overlapped_features features that overlap.
#' @return data.table
#' @keywords internal
.filterOverlapped = function(input, summary_function, overlapped_features) {
    overlapped = merge(input, overlapped_features, by = "feature")
    overlapped = overlapped[, .(agg_intensity = summary_function(Intensity, na.rm = TRUE)),
                            by = c("feature", "id")]
    maximized = overlapped[, .(max_intensity = max(agg_intensity)), by = "feature"]
    overlapped = merge(overlapped, maximized, by = "feature")
    overlapped = overlapped[overlapped$agg_intensity != overlapped$max_intensity] # drop == !
    # TODO: there is a better pattern for this
    input[!(input$id %in% overlapped$id), ]
}
