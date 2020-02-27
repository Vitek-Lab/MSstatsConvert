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

OpenMStoMSstatsFormat <- function(
  input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max) {
  
  fewMeasurements = .isLegalValue(fewMeasurements, c("remove", "keep"))
  input = .cleanRawOpenMS(input)
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
    input
  )
  
  feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon",
                   "ProductCharge")
  input = .handleSharedPeptides(input, useUniquePeptide)
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows, 
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature)
  input = merge(input[, setdiff(colnames(input), c("Condition", "BioReplicate")), 
                      with = FALSE], annotation, by = "Run")
  input # ProteinName as factor?
}
