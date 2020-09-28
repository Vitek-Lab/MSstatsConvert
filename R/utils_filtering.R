#' Handle PSM/proteins scores
#' @param input `data.table` preprocessed by one of the .cleanRaw* functions.
#' @param score_filtering list of by-score filtering controls.
#' @param exact_filtering list of exact filtering controls.
#' @param pattern_filtering list of by-pattern filtering controls.
#' @return data.table
#' @keywords internal
.handleFiltering = function(input, score_filtering, exact_filtering,
                            pattern_filtering) {
    for (filtering in score_filtering) {
        input = .filterByScore(input, filtering[["score_column"]], 
                               filtering[["score_threshold"]], 
                               filtering[["direction"]], filtering[["behavior"]], 
                               filtering[["handle_na"]], filtering[["fill_value"]],
                               filtering[["filter"]], filtering[["drop_column"]])
    }
    for (filtering in exact_filtering) {
        input = .filterExact(input, filtering[["col_name"]], 
                             filtering[["filter_symbols"]],
                             filtering[["behavior"]], filtering[["fill_value"]],
                             filtering[["filter"]], filtering[["drop_column"]])
    }
    for (filtering in pattern_filtering) {
        input = .filterByPattern(input, filtering[["col_name"]], 
                                 filtering[["pattern"]], 
                                 filtering[["filter"]], 
                                 filtering[["drop_column"]])
    }
    input
}


#' Filter PSMs / proteins by a given score column.
#' @param input `data.table` preprocessed by one of the .cleanRaw* functions.
#' @param score_column chr, name of the column that contains scores.
#' @param score_threshold num, values below or above this threshold will be
#' removed from the data.
#' @param direction chr, if "greater" only values above the threshold will be
#' retained, if "smaller" - below the threshold.
#' @param behavior chr, if "remove", values below/above the threshold will be 
#' removed, if "replace", they will be set to `fill_value`.
#' @param fill_value if `behavior` = "replace", values below/above the threshold
#' will be replaced with `fill_value`. Defaults to `NA`.
#' @param filter If TRUE, filtering will be performed.
#' @param drop if TRUE, `score_column` will be removed.
#' @return data.table
#' @keywords internal
.filterByScore = function(input, score_column, score_threshold, direction,
                          behavior, handle_na = "keep", fill_value = NA,
                          filter = TRUE, drop = TRUE) {
    if (!is.element(score_column, colnames(input))) {
        msg = paste("**", score_column, "not found in input columns.")
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        return(input)
    }
    
    if (filter) {
        if (direction == "greater") {
            score_filter = input[[score_column]] >= score_threshold
        } else {
            score_filter = input[[score_column]] <= score_threshold
        }
        
        if (handle_na == "keep") {
            score_filter = score_filter | is.na(input[[score_column]])
        } else{
            score_filter = score_filter & !is.na(input[[score_column]])
        }
        
        if (behavior == "remove") {
            input = input[score_filter, ]    
        } else {
            input[!score_filter, "Intensity"] = fill_value
        }
        
        msg = .makeScoreFilterMessage(score_column, score_threshold, direction,
                                      behavior, fill_value)
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    
    if (drop) {
        input = input[, !(colnames(input) == score_column), with = FALSE]
    }
    
    input
}


#' Handle filtering by pattern
#' @param input `data.table` preprocessed by one of the .cleanRaw* functions.
#' @param col_name chr, name of the column with peptide sequences.
#' @param pattern chr, regular expression - matching peptides will be 
#' removed from the data.
#' @param filter lgl, if TRUE, peptides will be actually filtered.
#' @param drop lgl, if TRUE, the `column` will be dropped.
#' @return data.table
#' @keywords internal
.filterByPattern = function(input, col_name, patterns, filter, drop) {
    if (!is.element(col_name, colnames(input))) {
        msg = paste(col_name, "not found in input columns.")
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        return(input)
    }
    
    if (filter) {
        msg = paste("** Sequences containing", 
                    paste(patterns, sep = ", ", collapse = ", "), "are removed.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        pattern_filter = rep(TRUE, nrow(input))
        for (pattern in patterns) {
            pattern_filter = pattern_filter & !grepl(pattern, input[[col_name]])
        }
    } else {
        pattern_filter = rep(TRUE, nrow(input))
    }
    
    if (drop) {
        input[pattern_filter, !(colnames(input) == col_name), with = FALSE]
    } else {
        input[pattern_filter, ]
    }
}


#' Filter out specified symbols.
#' @param input data.table preprocessed by one of the .cleanRaw* functions.
#' @param col_name chr, name of the column that will be the base for filtering
#' @param filter_symbols character vector of symbols that will be removed
#' @param behavior chr, if "remove", values below/above the threshold will be 
#' removed, if "replace", they will be set to `fill_value`.
#' @param fill_value if `behavior` = "replace", values below/above the threshold
#' will be replaced with `fill_value`. Defaults to `NA`.
#' @param filter lgl, if TRUE, decoy proteins will be removed from the data.
#' @param drop lgl, if TRUE, column that contains decoy proteins will be dropped.
#' @return data.table
#' @keywords internal
.filterExact = function(input, col_name, filter_symbols, behavior, 
                        fill_value, filter, drop
) {
    if (!is.element(col_name, colnames(input))) {
        msg = paste(col_name, "not found in input columns.")
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
        return(input)
    }
    
    find_col = colnames(input) == col_name
    if (filter) {
        exact_filter = !(input[[col_name]] %in% filter_symbols)
        msg = .makeExactFilterMessage(col_name, filter_symbols, 
                                      behavior, fill_value)
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    } else {
        exact_filter = rep(TRUE, nrow(input))
    }
    
    if (behavior == "remove") {
        input = input[exact_filter]
    } else {
        input$Intensity = ifelse(exact_filter, input$Intensity, fill_value)
    }
    
    
    if (drop) {
        input = input[, !find_col, with = FALSE]    
    }
    input
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
    for (col in intersect(filter_columns, colnames(input))) {
        has_col = is.element(col, colnames(input))
        input = .filterExact(input, col, filter_symbols, "remove", NULL,
                             has_col, has_col)        
    }
    input
}
