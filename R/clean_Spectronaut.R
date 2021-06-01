#' Clean raw Spectronaut output.
#' @param msstats_object an object of class `MSstatsSpectronautFiles`.
#' @param intensity chr, specifies which column will be used for Intensity.
#' @return `data.table`
#' @keywords internal
.cleanRawSpectronaut = function(msstats_object, intensity) {
  FFrgLossType = FExcludedFromQuantification = NULL
  
  spec_input = getInputFile(msstats_object, "input")
  spec_input = spec_input[FFrgLossType == "noloss", ]
  
  if (is.character(spec_input$FExcludedFromQuantification)) {
      spec_input = spec_input[FExcludedFromQuantification == "False", ]
  } else {
      spec_input = spec_input[!(as.logical(FExcludedFromQuantification)), ]
  }

  f_charge_col = .findAvailable(c("FCharge", "FFrgZ"), colnames(spec_input))
  pg_qval_col = .findAvailable(c("PGQvalue"), colnames(spec_input))
  cols = c("PGProteinGroups", "EGModifiedSequence", "FGCharge", "FFrgIon", 
           f_charge_col, "RFileName", "RCondition", "RReplicate", 
           "EGQvalue", pg_qval_col, paste0("F", intensity))
  cols = intersect(cols, colnames(spec_input))
  spec_input = spec_input[, cols, with = FALSE]
  data.table::setnames(
    spec_input, 
    c("PGProteinGroups", "EGModifiedSequence", "FGCharge", "FFrgIon",
      f_charge_col, "RFileName", paste0("F", intensity), "RCondition", "RReplicate"),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Intensity", "Condition", "BioReplicate"), 
    skip_absent = TRUE)
  .logSuccess("Spectronaut", "clean")
  spec_input
}
