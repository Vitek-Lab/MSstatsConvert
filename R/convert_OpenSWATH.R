#' Import OpenSWATH files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenSWATH, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param filter_with_mscore TRUE(default) will filter out the features that have greater than mscore_cutoff in m_score column. Those features will be removed.
#' @param mscore_cutoff Cutoff for m_score. Default is 0.01.
#'  
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 

OpenSWATHtoMSstatsFormat = function(
  input, annotation, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
  useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  .setMSstatsLogger(use_log_file, append, verbose)
  # .checkConverterParams()
  
  input = .cleanRawOpenSWATH(input)
  annotation = .makeAnnotation(input, .getDataTable(annotation))
  
  input = .handleFiltering(input, "m_score", mscore_cutoff, "smaller", "remove")
  input = .filterExact(input, "decoy", 1)
  input = .cleanRawOpenSWATH(input)
  input = .handleSharedPeptides(input, useUniquePeptide)
  feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon")
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature,
                                         feature_cols)
  input = .mergeAnnotation(input, annotation)
  input
}


#' Clean raw OpenSWATH files
#' @param os input OpenSWATH report or a path to it.
#' @param ... optional, additional parameters for data.table::fread
#' @return data.table
#' @keywords internal
.cleanRawOpenSWATH = function(os_input, ...) {
  os_input = .getDataTable(os_input, ...)
  colnames(os_input) = .updateColnames(
    os_input,
    c("FullPeptideName", "Charge", "filename"),
    c("PeptideSequence", "PrecursorCharge", "Run"))
  os_input = os_input[, lapply(.(aggr_Fragment_Annotation, aggr_Peak_Area), 
                               function(x) unlist(tstrsplit(x, ";", fixed = TRUE))),
                      by = c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run")]
  colnames(os_input) = .updateColnames(os_input, c("V1", "V2"), 
                                       c("FragmentIon", "Intensity"))
  os_input[, c("PeptideSequence", "FragmentIon")] = os_input[, lapply(.(PeptideSequence, FragmentIon), 
                                                                   function(x) gsub(":", "_", x))]
  os_input[["Intensity"]] = as.numeric(os_input[["Intensity"]])
  os_input[os_input[["Intensity"]] < 1, "Intensity"] = NA
  os_input = .fillValues(os_input, c("ProductCharge" = NA, "IsotopeLabelType" = "L"))
  os_input
}
