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

.OpenMSFromDFs = function(
  input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  
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


.OpenMSFromFiles = function(
  input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  
}

setGeneric("OpenMStoMSstatsFormat",
           function(input, annotation, ...) {
             standardGeneric("OpenMStoMSstatsFormat")
           })

setMethod("OpenMStoMSstatsFormat",
          signature = rep("character", 2),
          function(input, annotation, ...)
            .OpenMSFromFiles(input, annotation, ...))


setMethod("OpenMStoMSstatsFormat",
          signature = rep("data.frame", 2),
          function(input, annotation, ...) {
            .OpenMSFromDFs(input, annotation, ...)
          })

setMethod("OpenMStoMSstatsFormat",
          signature = c("data.frame", "missing"),
          function(input, ...) {
            .OpenMSFromDFs(input, ...)
          })

setMethod("OpenMStoMSstatsFormat",
          signature = c("character", "missing"),
          function(input, ...) {
            .OpenMSFromFiles(input, ...)
          })


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
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max
) {
  
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


.cleanRawOpenMS = function(om_input) {
  colnames(om_input) = .standardizeColnames(om_input)
  om_input = data.table::as.data.table(om_input)
  om_input[["Intensity"]] = as.numeric(om_input[["Intensity"]])
  if(!is.element("IsotopeLabelType", colnames(om_input))) {
    om_input = .fillValues(om_input, c("IsotopeLabelType" = "L"))
  }
  om_input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
               "FragmentIon", "ProductCharge", "IsotopeLabelType",
               "Condition", "BioReplicate", "Run", "Intensity"),
           with = FALSE]
}