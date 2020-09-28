#' Clean raw output from OpenMS
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @return data.table
#' @keywords internal
.cleanRawOpenMS = function(msstats_object) {
  Intensity = PSM = NULL
  
  om_input = getInputFile(msstats_object, "input")
  om_input[, Intensity := as.numeric(as.character(Intensity))]
  om_input = unique(om_input)
  
  if (getDataType(msstats_object) == "MSstats") {
    if (!is.element("IsotopeLabelType", colnames(om_input))) {
      om_input = .fillValues(om_input, c("IsotopeLabelType" = "L"))
    }
  } else {
    om_input[, PSM := do.call(".combine", .SD), 
             .SDcols = c("PeptideSequence", "Charge", "RetentionTime")]
    om_input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
  }
  
  all_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", "Charge",
               "FragmentIon", "ProductCharge", "PSM", "IsotopeLabelType",
               "Condition", "BioReplicate", "Run", "Channel", "Intensity",
               "Mixture", "TechRepMixture", "TechReplicate", "Fraction")
  cols = intersect(all_cols, colnames(om_input))
  om_input = om_input[, cols, with = FALSE]
  data.table::setnames(om_input, "Charge", "PrecursorCharge", 
                       skip_absent = TRUE)
  unique(om_input)
}