#' Remove peptides assigned to more than one protein.
#' @param input data.table pre-processed by one of the .cleanRaw* functions.
#' @param protein_column chr, name of the column with names of proteins.
#' @param peptide_column chr, name of the column with peptide sequences.
#' @return data.table
#' @keywords internal
.removeSharedPeptides = function(input, protein_column, peptide_column) {
    count = NULL
    
    unique_pairs = unique(input[, c(protein_column, peptide_column), 
                                with = FALSE])
    protein_counts = unique_pairs[, list(count = .N), by = peptide_column]
    protein_counts = unique(protein_counts[count == 1L, peptide_column, 
                                           with = FALSE])
    merge(input, protein_counts, sort = FALSE)[, colnames(input), with = FALSE]
}

#' Handle shared peptides.
#' @inheritParams .removeSharedPeptides
#' @param remove_shared lgl, if TRUE, shared peptides will be removed
#' @keywords internal
.handleSharedPeptides = function(input, remove_shared = TRUE,
                                 protein_column = "ProteinName",
                                 peptide_column = "PeptideSequence") {
    if (remove_shared) {
        input = .removeSharedPeptides(input, protein_column, peptide_column)
        getOption("MSstatsLog")("INFO", "Shared peptides are removed.")
        getOption("MSstatsMsg")("INFO", "Shared peptides are removed.")
    } 
    input
}
