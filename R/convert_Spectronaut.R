#' Import Spectronaut files
#' 
#' @param input name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Spectronaut, use annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area. 'NormalizedPeakArea' uses peak area normalized by Spectronaut.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 

SpectronauttoMSstatsFormat = function(
  input, annotation = NULL, intensity = 'PeakArea', filter_with_Qvalue = TRUE,
  qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements="remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  
  fewMeasurements = .isLegalValue(fewMeasurements, c("remove", "keep"))
  intensity = .isLegalValue(intensity, c("PeakArea", "NormalizedPeakArea"))
  .checkColumns("Input", 
                c("F.FrgLossType", "F.ExcludedFromQuantification",
                  "PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge",
                  "F.FrgIon", "R.FileName", "EG.Qvalue"), colnames(input))
  .checkColumns("Input", c("F.PeakArea", "F.NormalizedPeakArea"), 
                colnames(input), "optional")
  .checkColumns("Input", c("F.Charge", "F.FrgZ"), colnames(input), "optional")
  
  .setMSstatsLogger(use_log_file, append, verbose)
  input = .getDataTable(input)
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
    input, 
    c("R.FileName" = "Run", "R.Condition" = "Condition", "R.Replicate" = "BioReplicate")
  )  
  input = .cleanRawSpectronaut(input)
  input = .handleFiltering(input, "PG.Qvalue", 0.01, "greater", "fill", NA)
  # TODO: 1. Does 0.01 have to be hard-coded? 2. Explain in documentation that this is protein q-value. 3. Log+message
  input = .handleFiltering(input, "Qvalue", qvalue_cutoff, "greater", "fill", 0, TRUE)
  # TODO: 1. Explain in documentation that this is precursor q-value. 2. Log+message
  input = .handleSharedPeptides(input, useUniquePeptide)
  feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows, fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature)
  input = merge(input[, setdiff(colnames(input), 
                                c("Condition", "BioReplicate")), 
                      with = FALSE], annotation, by = "Run", all = TRUE)
  input = .fillValues(input, c("IsotopeLabelType" = "L"))
  input # Convert ProteinName to factor?
}


.cleanRawSpectronaut = function(spec_input, intensity) {
  colnames(spec_input) = .standardizeColnames(colnames(spec_input))
  
  spec_input = spec_input[spec_input[["F.FrgLossType"]] == "noloss", ]
  spec_input = spec_input[!spec_input[["F.ExcludedFromQuantification"]], ]
  # XIC quality. TODO: explain in documentation
  
  f_charge_col = .findAvailable(c("F.Charge", "F.FrgZ"), colnames(spec_input))
  pg_qval_col = .findAvailable(c("PG.Qvalue"), colnames(spec_input))
  cols = c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon", 
           f_charge_col, "R.FileName", "EG.Qvalue", pg_qval_col, 
           paste0("F.", intensity))
  spec_input = spec_input[, cols, with = FALSE]
  colnames(spec_input) = .updateColnames(
    spec_input, 
    c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon",
      f_charge_col, "R.FileName", "EG.Qvalue", paste0("F.", intensity)),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Qvalue", "Intensity"))
  spec_input
}
