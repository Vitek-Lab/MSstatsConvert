#' Import Skyline files
#'
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 

SkylinetoMSstatsFormat = function(
    input, annotation = NULL, removeiRT = TRUE, filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE,
    use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
    .setMSstatsLogger(use_log_file, append, verbose)
    # .checkConverterParams()
    
    input = .cleanRawSkyline(input)
    annotation = .makeAnnotation(input, annotation)
    
    input = .filterExact(input, "ProteinName", 
                         c("DECOY", "Decoys"), FALSE)
    input = .filterExact(input, "StandardType", "iRT", FALSE, removeiRT)
    input = .handleOxidationPeptides(input, "PeptideSequence", "+16", 
                                     removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .handleFiltering(input, "Truncated", 0L, "smaller", "fill", NA_real_)
    input = .handleIsotopicPeaks(input)
    input = .handleFiltering(input, "DetectionQValue", qvalue_cutoff, "smaller", 
                             "fill", 0, TRUE, filter_with_Qvalue)
    # TODO: conditionally check if DetectionQValue is in the input at the beginning
    feature_cols = c("PeptideSequence", "PrecursorCharge")
    input = .cleanByFeature(input, feature_cols, max, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature,
                                           feature_cols)
    input = .mergeAnnotation(input, annotation)
    input = .fillValues(input, c("ProductCharge" = NA, "IsotopeLabelType" = "L",
                                 "FragmentIon" = "sum"))
    .MSstatsFormat(input)
}


#' Clean raw data from Skyline
#' @param sl_input data.frame with raw Skyline input or a path to Skyline output.
#' @param ... optional, additional parameters to data.table::fread.
#' @return data.table
#' @keywords internal
.cleanRawSkyline = function(sl_input, ...) {
    sl_input = .getDataTable(sl_input, ...)
    colnames(sl_input) = gsub("\\.", "", colnames(sl_input))
    colnames(sl_input) = .updateColnames(sl_input, c("FileName", "Area"),
                                         c("Run", "Intensity"))
    colnames(sl_input) = .standardizeColnames(colnames(sl_input))
    
    sl_input = sl_input[, !(colnames(sl_input) == "PeptideSequence"), with = FALSE]
    colnames(sl_input) = .updateColnames(sl_input, "PeptideModifiedSequence",
                                         "PeptideSequence")
    sl_input[["Intensity"]] = as.numeric(sl_input[["Intensity"]])
    if (is.element("DetectionQValue", colnames(sl_input))) {
        sl_input[["DetectionQValue"]] = as.numeric(as.character(sl_input[["DetectionQValue"]]))    
    }
    if(is.character(sl_input[["Truncated"]])) {
        sl_input[["Truncated"]] = sl_input[["Truncated"]] == "True"
    }
    sl_input[["Truncated"]] = as.integer(sl_input[["Truncated"]])
    
    sl_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition",
                "BioReplicate", "Run", "Intensity", "StandardType")
    sl_cols = c(sl_cols, "Fraction", "DetectionQValue", "Truncated")
    sl_input = sl_input[, intersect(sl_cols, colnames(sl_input)), with = FALSE]
    non_missing = sapply(sl_input, function(x) !all(is.na(x))) | colnames(sl_input) == "StandardType"
    sl_input[, non_missing, with = FALSE]
}


#' Handle isotopic peaks
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @return data.table
#' @keywords internal
.handleIsotopicPeaks = function(input) {
    if (.checkDDA(input)) {
        input = .summarizeMultipleMeasurements(input, 
                                               function(x) sum(x, na.rm = TRUE),
                                               c("ProteinName",
                                                 "PeptideSequence",
                                                 "PrecursorCharge",
                                                 "Run")) # or just PeptideSequence and PrecursorCharge?
    }
    getOption("MSstatsLog")("INFO", "Three isotopic preaks per feature and run are summed")
    getOption("MSstatsMsg")("INFO", "Three isotopic preaks per feature and run are summed")
    input
}


#' Check validity of DDA data
#' @inherit .handleIsotopicPeaks
#' @return logical
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
        msg = paste("Please check precursors information.",
                    "If your experiment is DIA, please remove the precursors.",
                    "If your experiments is DDA, please check the precursor information.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    length(check_DDA) < 3
}
