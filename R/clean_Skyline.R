#' Clean raw data from Skyline
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @return data.table
#' @keywords internal
.cleanRawSkyline = function(msstats_object) {
    Truncated = isZero = DetectionQValue = Intensity = NULL
    
    sl_input = getInputFile(msstats_object, "input")
    data.table::setnames(sl_input, c("FileName", "Area"),
                         c("Run", "Intensity"), skip_absent = TRUE)
    
    sl_input = sl_input[, !(colnames(sl_input) == "PeptideSequence"), 
                        with = FALSE]
    data.table::setnames(sl_input, "PeptideModifiedSequence",
                         "PeptideSequence", skip_absent = TRUE)
    suppressWarnings({
        sl_input[, Intensity := as.numeric(as.character(Intensity))]
    })
    if (is.element("DetectionQValue", colnames(sl_input))) {
        sl_input[, DetectionQValue := as.numeric(as.character(DetectionQValue))]
    }
    if (is.character(sl_input$Truncated) | is.factor(sl_input$Truncated)) {
        sl_input[, Truncated := as.character(Truncated) == "True"]
    }
    sl_input[, Truncated := as.character(Truncated)]
    
    sl_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition",
                "BioReplicate", "Run", "Intensity", "StandardType")
    sl_cols = c(sl_cols, "Fraction", "DetectionQValue", "Truncated")
    sl_input = sl_input[, intersect(sl_cols, colnames(sl_input)), with = FALSE]
    sl_input[, isZero := ((Intensity == 0) | (Intensity < 1)) & !is.na(Intensity)]
    .logSuccess("Skyline", "clean")
    sl_input
}


#' Handle isotopic peaks
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @param aggregate if TRUE, isotopic peaks will be summed.
#' @return data.table
#' @keywords internal
.handleIsotopicPeaks = function(input, aggregate = FALSE) {
    FragmentIon = ProductCharge = NULL
    
    if (aggregate) {
        input$isZero = ifelse(input$Intensity == 0, TRUE, input$isZero)
        if (.checkDDA(input)) {
            feature_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                             "Run")
            feature_cols = c(feature_cols, c("BioReplicate", "Condition", "StandardType", 
                                             "IsotopeLabelType", "DetectionQValue"))
            feature_cols = intersect(feature_cols, colnames(input))
            input = .summarizeMultipleMeasurements(input, sum, feature_cols)
            input[, FragmentIon := NA]
            input[, ProductCharge := NA]
            getOption("MSstatsLog")("INFO", "** Three isotopic preaks per feature and run are summed")
            getOption("MSstatsMsg")("INFO", "** Three isotopic preaks per feature and run are summed")
        }
    }
    input
}


#' Check validity of DDA data
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @return logical
#' @return logical, TRUE means that the `input` dataset comes from a DDA experiment
#' @keywords internal
.checkDDA = function(input) {
    fragment_ions = as.character(unique(input[["FragmentIon"]]))
    check_DDA = setdiff(c("precursor", "precursor [M+1]", "precursor [M+2]"), 
                        fragment_ions)
    frags = setdiff(fragment_ions, 
                    c('precursor', 'precursor [M+1]', 'precursor [M+2]'))
    precursors = intersect(fragment_ions, 
                           c("precursor", "precursor [M+1]", "precursor [M+2]"))
    if (length(frags) > 0 & length(precursors) > 0) {
        msg = paste("** Please check precursors information in the FragmentIon column.",
                    "If your experiment is DIA or SRM, please remove the precursors.",
                    "If your experiments is DDA, please check the precursor information.",
                    "MSstats only performs analysis on either MS1 or MS2, but not both simultaneously.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    length(check_DDA) < 3
}
