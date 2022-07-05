#' Clean raw DIA-NN output
#' @param msstats_object Object that inherits from MSstatsInputFiles class.
#' @param protein_column chr, name of a column that will be used as protein ID.
#' @param fragment_column chr, name of a column that will be used as a source
#' of fragment intensities.
#' @return data.table
#' @keywords internal
.cleanRawDIANN = function(msstats_object, protein_column,
                          fragment_column) {
  diann_input = getInputFile(msstats_object, "input")
  diann_input = diann_input[, c("Run", protein_column, 
                                "Modified.Sequence", "Precursor.Charge",
                                "Q.Value", "Protein.Q.Value", 
                                fragment_column), with = FALSE]
  setnames(diann_input, c("Run", "ProteinName", "PeptideSequence",
                          "PrecursorCharge", "Qvalue", "ProteinQValue",
                          "FragmentQuants"))
  diann_fragments = diann_input[, tstrsplit(FragmentQuants, ";")]
  setnames(diann_fragments, paste0("Frag", seq_len(ncol(diann_fragments))))
  diann_input[, FragmentQuants := NULL]
  diann_input = cbind(diann_input, diann_fragments)
  fragment_quant_cols = colnames(diann_input)[grepl("Frag", 
                                                    colnames(diann_input))]
  diann_input = melt(diann_input, measure.vars = fragment_quant_cols)
  setnames(diann_input, c("variable", "value"), c("FragmentIon", "Intensity"))
  diann_input = diann_input[!is.na(Intensity)]
  diann_input[, ProductCharge := 1]  
  diann_input
}
