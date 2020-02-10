#' A dummy function to store shared documentation items.
#' 
#' @param fewMeasurements 'remove'(default) will remove the features that have 1 or 2 measurements across runs.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.
#' @param removeMpeptides TRUE will remove the peptides including 'M' sequence. FALSE is default.
#' 
#' @keywords internal
#' 

.documentFunction = function(fewMeasurements, 
                             useUniquePeptide,
                             summaryforMultipleRows, 
                             removeProtein_with1Feature,
                             removeProtein_with1Protein,
                             removeOxidationMpeptides,
                             removeMpeptides) {
    
}

.updateColnames = function(data_frame, column_update) {
    columns <- colnames(data_frame)
    not_changing <- setdiff(columns, names(column_update))
    column_update[not_changing] <- not_changing
    unname(column_update[columns])
}


# ALL DATA CHECKS!!

.isLegalValue = function(parameter, legal_values = NULL, 
                         can_be_null = FALSE) {
    parameter_name = deparse(substitute(parameter))
    if(is.null(parameter)) {
        if(!can_be_null) {
            stop(paste("Parameter", parameter_name, "cannot be NULL"))    
        }
    } else {
        if(!is.null(legal_values)) {
            if(!is.element(parameter, legal_values)) {
                stop(paste("Parameter", parameter_name, "must be one of", 
                           paste(legal_values, sep = ", ", collapse = ", ")))
            }    
        }
    }
    invisible(TRUE)
}

# COMMON DATA PREPROCESSING

.removeSharedPeptides = function(data_frame, proteins_column, peptides_column) {
    unique_pairs = unique(data_frame[, c(proteins_column, peptides_column)])
    protein_counts = aggregate(x = unique_pairs[[proteins_column]], 
                               by = list(peptide = unique_pairs[[peptides_column]]),
                               length)
    counts = protein_counts[["x"]]
    names(counts) = protein_counts[["peptide"]]
    if(length(counts) == 0) {
        data_frame
    } else {
        data_frame[counts[data_frame[[peptides_column]]] == 1, ]    
    }
    # TODO: message for the user / log
}
