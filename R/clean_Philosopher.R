#' Clean raw Philosopher files
#' @param msstats_object object of class MSstatsPhilosopherFiles
#' @param protein_id_col character name of a column that identifies proteins
#' @param peptide_id_col character name of a column that identifies peptides
#' @param channels character vector of channel labels
#' @param remove_shared_peptides logical, if TRUE, shared peptides will be 
#' removed based on the IsUnique column from Philosopher output
#' @keywords internal
#' @return data.table
.cleanRawPhilosopher = function(msstats_object, protein_id_col,
                                peptide_id_col, channels,
                                remove_shared_peptides) {
  ProteinName = ModifiedPeptideSequence = PeptideSequence = NULL
  PrecursorCharge = Intensity = IsUnique = PSM = NULL
  
  channels = .standardizeColnames(channels)
  input = getInputFile(msstats_object, "psm_data")  
  protein_col = .standardizeColnames(protein_id_col)
  peptide_col = .standardizeColnames(peptide_id_col)
  
  input[, ModifiedPeptideSequence := ifelse(is.na(ModifiedPeptideSequence),
                                            PeptideSequence, 
                                            ModifiedPeptideSequence)]
  data.table::setnames(input, stringi::stri_replace_all(colnames(input), 
                                                        fixed = "Channel", replacement = ""))
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
  input[, PSM := paste(PeptideSequence, PrecursorCharge, 
                       1:nrow(input), sep = "_")]
  input = data.table::melt(input, measure.vars = channels,
                           variable.name = "Channel", value.name = "Intensity",
                           variable.factor = FALSE)
  input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
  input[, Run := stringi::stri_split(Run, regex = "\\.", simplify = TRUE)[, 1]]
  input = input[!is.na(ProteinName) & ProteinName != ""]
  input
}
