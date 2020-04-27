#' Import OpenMS files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenMS, which includes feature(peptide ion)-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 

OpenMStoMSstatsFormat = function(
  input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
  .setMSstatsLogger(use_log_file, append, verbose)
  # .checkConverterParams()

  input = .cleanRawOpenMS(input)
  annotation = .makeAnnotation(input, .getDataTable(annotation))

  feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon",
                   "ProductCharge")
  input = .handleSharedPeptides(input, useUniquePeptide)
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows, 
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature,
                                         feature_cols)
  input = .mergeAnnotation(input, annotation)
  input
}


#' Clean raw output from OpenMS
#' @param om_input OpenMS report or a path to it.
#' @param ... optional, additional parameters for data.table::fread.
#' @return data.table
#' @keywords internal
.cleanRawOpenMS = function(om_input, ...) {
  om_input = .getDataTable(om_input, ...)
  colnames(om_input) = .standardizeColnames(colnames(om_input))
  om_input[["Intensity"]] = as.numeric(om_input[["Intensity"]])
  if (!is.element("IsotopeLabelType", colnames(om_input))) {
    om_input = .fillValues(om_input, c("IsotopeLabelType" = "L"))
  }
  om_input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
               "FragmentIon", "ProductCharge", "IsotopeLabelType",
               "Condition", "BioReplicate", "Run", "Intensity"),
           with = FALSE]
}


# TODO: modify parameter names 
OpenMStoMSstatsTMTFormat = function(
  input, useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
  rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
  summaryforMultiplePSMs = sum, use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
  .setMSstatsLogger(use_log_file, append, verbose)
  # .checkConverterParams()
  
  input = .cleanRawOpenMSTMT(input)
  feature_cols = c("PeptideSequence", "Charge", "Reference", "RetentionTime") # THIS used to include protein name and run
  input = .removeMissingAllChannels(input, feature_cols)
  input = .handleSharedPeptides(input, useUniquePeptide)
  input = .cleanByFeatureTMT(input, feature_cols, summaryforMultiplePSMs, 
                             rmPSM_withfewMea_withinRun, rmPSM_withMissing_withinRun)
  input = .handleSingleFeaturePerProtein(input, rmProtein_with1Feature,
                                         c("PeptideSequence", "Charge")) # or just PSM
  input = .handleFractions(input)
  input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
            "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Intensity")]
}


#' Clean raw OpenMS TMT data
#' @param om_input OpenMS report or a path to it.
#' @param ... optional, additional parameters for data.table::fread
#' @return data.table
#' @keywords internal
.cleanRawOpenMSTMT = function(om_input, ...) {
  om_input = .getDataTable(om_input, ...)
  colnames(om_input) = .standardizeColnames(colnames(om_input))
  om_input$PSM = paste(om_input$PSM, 1:nrow(om_input), sep = "_")
  om_input
}
