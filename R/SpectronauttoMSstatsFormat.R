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
SpectronauttoMSstatsFormat = function(
  input, annotation = NULL, intensity = 'PeakArea', filter_with_Qvalue = TRUE,
  qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements="remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max) {
  
  .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))
  .isLegalValue(intensity, legal_values = c("PeakArea", "NormalizedPeakArea"))
  .checkColumns("Input", 
                c("F.FrgLossType", "F.ExcludedFromQuantification",
                  "PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge",
                  "F.FrgIon", "R.FileName", "EG.Qvalue"), colnames(input))
  .checkColumns("Input", c("F.PeakArea", "F.NormalizedPeakArea"), 
                colnames(input), "optional")
  .checkColumns("Input", c("F.Charge", "F.FrgZ"), colnames(input), "optional")
  
  input = input[input[["F.FrgLossType"]] == "noloss", ]
  input = input[!input[["F.ExcludedFromQuantification"]], ] # XIC quality. TODO: explain in documentation
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
    input, 
    c("R.FileName" = "Run", "R.Condition" = "Condition", "R.Replicate" = "BioReplicate")
  )  
  
  f_charge_col = .findAvailable(c("F.Charge", "F.FrgZ"), colnames(input))
  pg_qval_col = .findAvailable(c("PG.Qvalue"), colnames(input))
  input = .selectColumns(
    input, 
    c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon", 
      f_charge_col, "R.FileName", "EG.Qvalue", pg_qval_col, paste0("F.", intensity)))
  colnames(input) = .updateColnames(
    input, c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon",
             f_charge_col, "R.FileName", "EG.Qvalue", paste0("F.", intensity)),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Qvalue", "Intensity"))
  
  input = .handleFiltering(input, "PG.Qvalue", 0.01, "greater", "fill", NA)
  # TODO: 1. Does 0.01 have to be hard-coded? 2. Explain in documentation that this is protein q-value. 3. Log+message
  input = .handleFiltering(input, "Qvalue", qvalue_cutoff, "greater", "fill", 0, TRUE)
  # TODO: 1. Explain in documentation that this is precursor q-value. 2. Log+message
  input = .handleSharedPeptides(input, "ProteinName", "PeptideSequence",
                                remove_shared = useUniquePeptide)
  input = .cleanByFeature(input, c("PeptideSequence", "PrecursorCharge",
                                   "FragmentIon", "ProductCharge"),
                          summaryforMultipleRows, fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, 
                                         c("PeptideSequence", "PrecursorCharge",
                                           "FragmentIon", "ProductCharge"),
                                         removeProtein_with1Feature)
  # Conditional filtering of single-feature proteins.
  input <- merge(input, annotation, by = "Run", all = TRUE)
  input = .selectColumns(
    input, c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
             "ProductCharge", "Condition", "BioReplicate", "Run", "Intensity"))
  input = .fillValues(input, c("IsotopeLabelType" = "L"))
  .fixColumnTypes(input, factor_columns = "ProteinName")
}
