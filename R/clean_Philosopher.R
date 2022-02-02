.cleanRawPhilosopher = function(msstats_object, protein_id_col,
                                peptide_id_col, channels,
                                remove_shared_peptides) {
  channels = MSstatsConvert:::.standardizeColnames(channels)
  input = getInputFile(msstats_object, "psm_data")  
  protein_col = MSstatsConvert:::.standardizeColnames(protein_id_col)
  peptide_col = MSstatsConvert:::.standardizeColnames(peptide_id_col)
  
  input[, ModifiedPeptideSequence := ifelse(is.na(ModifiedPeptideSequence),
                                            PeptideSequence, 
                                            ModifiedPeptideSequence)]
  cols = c(protein_col, peptide_col, "Charge", "SpectrumName", 
           "PeptideProphetProbability", "Purity", "Modifications",
           "IsUnique", channels)
  cols = intersect(cols, colnames(input))
  input = input[, cols, with = FALSE]
  data.table::setnames(input, c(protein_col, peptide_col, "SpectrumName", "Charge"),
                       c("ProteinName", "PeptideSequence", "Run", "PrecursorCharge"))
  if (remove_shared_peptides) {
    input = input[(IsUnique)]
    input[, IsUnique := NULL]
  }
  input = data.table::melt(input, measure.vars = channels,
                           variable.name = "Channel", value.name = "Intensity",
                           variable.factor = FALSE)
  input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
  input = input[!is.na(ProteinName) & ProteinName != ""]
  input
}
