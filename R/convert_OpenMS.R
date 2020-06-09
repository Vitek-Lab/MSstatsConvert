#' Clean raw output from OpenMS
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @return data.table
#' @keywords internal
.cleanRawOpenMS = function(msstats_object) {
  om_input = getInputFile(msstats_object, "input")
  om_input[["Intensity"]] = as.numeric(om_input[["Intensity"]])
  
  if (getDataType(msstats_object) == "MSstats") {
    if (!is.element("IsotopeLabelType", colnames(om_input))) {
      om_input = .fillValues(om_input, c("IsotopeLabelType" = "L"))
    }
  } else {
    om_input$PSM = paste(om_input$PSM, 1:nrow(om_input), sep = "_")  
  }
  
  all_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", "Charge",
               "FragmentIon", "ProductCharge", "PSM", "IsotopeLabelType",
               "Condition", "BioReplicate", "Run", "Channel", "Intensity",
               "Fraction", "Reference", "RetentionTime")
  cols = intersect(all_cols, colnames(om_input))
  om_input = om_input[, cols, with = FALSE]
  colnames(om_input) = .updateColnames(om_input, "Charge", "PrecursorCharge")
  om_input
}