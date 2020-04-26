# TODO: modify parameter names
SpectroMinetoMSstatsTMTFormat <- function(
  input, annotation, filter_with_Qvalue = TRUE, qvalue_cutoff = 0.01,
  useUniquePeptide = TRUE, rmPSM_withMissing_withinRun = FALSE,
  rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE,
  summaryforMultipleRows = sum, use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
  
  .setMSstatsLogger(use_log_file, append, verbose)
  # annotation = .makeAnnotation() ...
  # TODO: checks go here
  input = .cleanRawSpectroMineTMT(input)
  input = .handleFiltering(input, "PQ.QValue", 0.01, "smaller", "fill", NA,
                           TRUE,) 
  input = .handleFiltering(input, "Qvalue", qvalue_cutoff, "smaller", "fill", 
                           NA, TRUE, filter_with_Qvalue)
  feature_cols = c("PeptideSequence", "Charge")
  input = .removeMissingAllChannels(input, feature_cols)
  input = .handleSharedPeptides(input, useUniquePeptide)
  input = .cleanByFeatureTMT(input, feature_cols, summaryforMultipleRows,
                             rmPSM_withfewMea_withinRun, rmPSM_withMissing_withinRun)
  
  input = .mergeAnnotation(input, annotation)
  input = .handleSingleFeaturePerProtein(input, rmProtein_with1Feature, "PSM")
  input = .handleFractions(input, annotation)
  input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
            "TechRepMixture", "Run", "Channel", "BioReplicate", "Condition", "Intensity")]
}

#' Clean raw SpectroMine TMT data
#' @param sm_input SpectroMine report or a path to it.
#' @param ... optional, additional parameters for data.table::fread
#' @importFrom data.table melt
#' @return data.table
#' @keywords internal
.cleanRawSpectroMineTMT = function(sm_input, ...) {
  sm_input = .getDataTable(sm_input, ...)
  colnames(sm_input) = .standardizeColnames(colnames(sm_input))
  # TODO: allow for other regular expressions?
  channels = .getChannelColumns(colnames(sm_input), "PSM", "Raw")
  if (length(channels) == 0L) {
    msg = paste("There is no channel intensity column in the input data,", 
                "which should start with 'PSM' and end with 'Raw'.")
    getOption("MSstatsMsg")("ERROR", msg)
    stop(msg)
  }
  sm_input = sm_input[, c("PG.ProteinAccessions", "P.MoleculeID", "PP.Charge",
                          "PG.QValue", "PSM.Qvalue", "R.FileName", channels),
                      with = FALSE] # TODO: check column validity
  colnames(sm_input) = .updateColnames(sm_input, 
                                       c("PG.ProteinAccessions", "P.MoleculeID", 
                                         "PP.Charge", "PSM.Qvalue", "R.FileName"),
                                       c("ProteinName", "PeptideSequence", "Charge",
                                         "Qvalue", "Run"))
  sm_input = sm_input[(sm_input$ProteinName != "") & (!is.na(sm_input$ProteinName)), ]
  sm_input = melt(sm_input, measure.vars = channels,
                  id.vars = setdiff(colnames(sm_input), channels),
                  variable.name = "Channel", value.name = "Intensity")
  sm_input$Channel = gsub("PSM", "", sm_input$Channel)
  sm_input$Channel = gsub("Raw", "", sm_input$Channel)
  sm_input$Channel = gsub(".", "", sm_input$Channel, fixed = TRUE)
  sm_input$PSM = paste(sm_input$PeptideSequence, sm_input$Charge, sep = "_")
  sm_input
}
