#' Import Progenesis files
#' 
#' @inheritParams .documentFunction
#' @param input name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.
#'
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek, Ulrich Omasits
#' 
#' @export
#' 

ProgenesistoMSstatsFormat <- function(
  input, annotation, useUniquePeptide = TRUE, summaryforMultipleRows = max,
  fewMeasurements = "remove", removeOxidationMpeptides = FALSE, 
  removeProtein_with1Peptide = FALSE) {
  fewMeasurements = .isLegalValue(fewMeasurements, c("remove", "keep"))
  # Check go here
  annotation = .makeAnnotation(annotation, 
                               c("Run" = "Run", "Condition" = "Condition",
                                 "BioReplicate" = "BioReplicate"))
  
  input = .cleanRawProgenesis(input, unique(annotation[["Run"]]))
  
  feature_cols = c("PeptideModifiedSequence", "PrecursorCharge")
  input = .handleOxidationPeptides(input, "PeptideModifiedSequence", 
                                   "Oxidation", removeOxidationMpeptides)
  input = .handleSharedPeptides(input, useUniquePeptide,
                                peptide_column = "PeptideModifiedSequence")
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, feature_cols, 
                                         removeProtein_with1Peptide)
  input = merge(input, annotation, by = "Run")
  input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
  input
}

