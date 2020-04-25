#' Filter PSMs / proteins by a given score column.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param score_column chr, name of the column that contains scores.
#' @param score_threshold num, values below or above this threshold will be
#' removed from the data.
#' @param direction chr, if "greater" only values above the threshold will be
#' retained, if "smaller" - below the threshold.
#' @param behavior chr, if "remove", values below/above the threshold will be 
#' removed, if "replace", they will be set to `fill_value`.
#' @param fill_value if `behavior` = "replace", values below/above the threshold
#' will be replaced with `fill_value`. Defaults to `NA`.
.filterByScore = function(input, score_column, score_threshold, direction,
                          behavior, fill_value = NA) {
    if (direction == "greater") {
        score_filter = input[[score_column]] >= score_threshold
    } else {
        score_filter = input[[score_column]] <= score_threshold
    }
    score_filter = score_filter & !is.na(input[[score_column]])
    if (behavior == "remove") {
        input = input[score_filter, ]    
    } else {
        input[!score_filter, "Intensity"] = fill_value
    }
    input
}


#' Handle PSM/proteins scores
#' @inheritParams .filterByScore
#' @param drop_column lgl, if TRUE, score column will be removed from the data.
#' @param filter lgl, if TRUE, data will be actually filtered.
#' @return data.table
#' @keywords internal
.handleFiltering = function(input, score_column, score_threshold, 
                            direction, behavior, fill_value = NA, 
                            drop_column = TRUE, filter = TRUE) {
    if (filter) { 
        input = .filterByScore(input, score_column, score_threshold, 
                               direction, behavior, fill_value)
    }
    if (drop_column) {
        input = input[, colnames(input) != score_column, with = FALSE] 
    }
    input
}


#' Handle oxidation or M-peptides.
#' @param data.table preprocessed by one of the .cleanRaw* functions.
#' @param sequence_column chr, name of the column with peptide sequences.
#' @param oxidation_pattern chr, regular expression - matching peptides will be 
#' removed from the data.
#' @param remove lgl, if TRUE, peptides will be actually filtered.
#' @return data.table
#' @keywords internal
.handleOxidationPeptides = function(input, sequence_column, 
                                    oxidation_pattern, remove) {
    if (remove) {
        msg = paste("Peptides containing", oxidation_pattern, "are removed")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        input[!(grepl(oxidation_pattern, input[[sequence_column]])), ]   
    } else {
        input
    }
}


#' Handle decoy protein.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param decoy_column chr, name of the column that contains names of decoy 
#' proteins.
#' @param drop lgl, if TRUE, column that contains decoy proteins will be dropped.
#' @param filter lgl, if TRUE, decoy proteins will be removed from the data.
#' @return data.table
#' @keywords internal
.filterExact = function(input, col_name, filter_symbols, 
                        drop = TRUE, filter = TRUE) {
    find_col = colnames(input) == col_name
    if (filter) {
        exact_filter = !(input[[col_name]] %in% filter_symbols)
    } else {
        exact_filter = rep(TRUE, nrow(input))
    }
    
    if (drop) {
        input[exact_filter, !find_col, with = FALSE]    
    } else {
        input[exact_filter, ]
    }
}


#' Filter rows that contain specifed symbols in multiple columns.
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @param filter_columns chr, names of columns in which elements will be matched
#' and removed.
#' @param filter_symbols chr, vector of strings. Rows with corresponding elements 
#' in `filter_columns` will be removed.
#' @return data.table
#' @keywords internal
.filterManyColumns = function(input, filter_columns, filter_symbols) {
    for(col in filter_columns) {
        has_col = is.element(col, colnames(input))
        input = .handleDecoyProteins(input, col, filter_symbols, 
                                     has_col, has_col)        
    }
    input
}
