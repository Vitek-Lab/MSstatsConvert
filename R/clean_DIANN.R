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
  cols = c("Run", protein_column, 
           "Modified.Sequence", "Precursor.Charge",
           "Q.Value", "Protein.Q.Value",  "Fragment.Info",
           fragment_column)
  
  diann_input = diann
  cols = intersect(cols, colnames(diann_input))
  diann_input = diann_input[, cols, with = FALSE]
  
  
  data.table::setnames(diann_input, cols,
                       MSstatsConvert:::.standardizeColnames)
  protein_column = MSstatsConvert:::.standardizeColnames(protein_column)
  fragment_column = MSstatsConvert:::.standardizeColnames(fragment_column)
  setnames(diann_input, 
           c(protein_column, "ModifiedSequence", "QValue",
             fragment_column),
           c("ProteinName", "PeptideSequence", "Qvalue",
             "FragmentQuants"))
  
  if (!is.element("FragmentInfo", colnames(diann_input))) {
    diann_fragments = diann_input[, tstrsplit(FragmentQuants, ";")]
    setnames(diann_fragments, paste0("Frag", seq_len(ncol(diann_fragments))))
    diann_input[, FragmentQuants := NULL]
    diann_input = cbind(diann_input, diann_fragments)
    fragment_quant_cols = colnames(diann_input)[grepl("Frag", 
                                                      colnames(diann_input))]
    diann_input = melt(diann_input, measure.vars = fragment_quant_cols)
    setnames(diann_input, c("variable", "value"), c("FragmentIon", "Intensity"))
  } else {
    diann_fragments = diann_input[, tstrsplit(FragmentQuants, ";")]
    diann_fragment_seqs = diann_input[, tstrsplit(FragmentInfo, ";")]
    diann_input[, FragmentQuants := NULL]
    diann_input[, FragmentInfo := NULL]
    data.table::setnames(diann_fragments, paste0("Frag", seq_len(ncol(diann_fragments))))
    data.table::setnames(diann_fragment_seqs, paste0("FragSeq", seq_len(ncol(diann_fragment_seqs))))
    diann_input = cbind(diann_input, diann_fragments, diann_fragment_seqs)
    
    diann_input = melt(diann_input, 
                       measure = list(colnames(diann_fragments),
                                      colnames(diann_fragment_seqs)),
                       value.name = c("Intensity", "FragmentIon"))
    diann_input[, variable := NULL]
    diann_input = diann_input[!is.na(FragmentIon) & !is.na(Intensity)]
  }
  diann_input = diann_input[!is.na(Intensity)]
  diann_input[, ProductCharge := 1]  
  diann_input
}
