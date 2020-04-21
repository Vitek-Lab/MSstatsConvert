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

.OpenSWATHFromDFs = function(
  input, annotation, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
  useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  
  # TODO: 1. Improved checks 2. Messages + logs
  fewMeasurements = .isLegalValue(fewMeasurements, 
                                  legal_values = c("remove", "keep"))
  annotation = .makeAnnotation(annotation, 
                               c("Run" = "Run", "Condition" = "Condition", 
                                 "BioReplicate" = "BioReplicate"))
  feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon")
  
  input = .handleFiltering(input, "m_score", mscore_cutoff, "smaller", "remove")
  input = .handleDecoyProteins(input, "decoy", 1)
  input = .cleanRawOpenSWATH(input)
  input = .handleSharedPeptides(input, useUniquePeptide)
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature)
  input = merge(input, annotation, by = "Run")
  input
}


.cleanRawOpenSWATH = function(openswath_input) {
  os_cols = c("ProteinName", "FullPeptideName", "Charge", "filename", 
              "aggr_Fragment_Annotation", "aggr_Peak_Area")
  input = data.table::as.data.table(openswath_input[, os_cols])
  colnames(input) = .updateColnames(
    input,
    c("FullPeptideName", "Charge", "filename"),
    c("PeptideSequence", "PrecursorCharge", "Run"))
  input = input[, lapply(.(aggr_Fragment_Annotation, aggr_Peak_Area), 
                         function(x) unlist(tstrsplit(x, ";", fixed = TRUE))),
                by = .(ProteinName, PeptideSequence, PrecursorCharge, Run)]
  colnames(input) = .updateColnames(input, c("V1", "V2"), 
                                    c("FragmentIon", "Intensity"))
  input[, c("PeptideSequence", "FragmentIon")] = input[, lapply(.(PeptideSequence, FragmentIon), 
                                                                function(x) gsub(":", "_", x))]
  input[["Intensity"]] = as.numeric(input[["Intensity"]])
  input[input[["Intensity"]] < 1, "Intensity"] = NA
  input = .fillValues(input, c("ProductCharge" = NA, "IsotopeLabelType" = "L"))
  input
}