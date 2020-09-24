#' Read file from a provided path or convert given data.frame to data.table
#' @param input report from a signal processing tool or a path to it
#' @param ... additional parameters for data.table::fread
#' @importFrom data.table as.data.table fread setDT
#' @keywords internal
.getDataTable = function(input, ...) {
    checkmate::checkTRUE(is.character(input) | inherits(input, "data.frame"))
    if (inherits(input, "data.frame")) {
        input = data.table::setDT(input)
    } else {
        input = data.table::fread(input, showProgress = FALSE, ...)
    }
    colnames(input) = .standardizeColnames(colnames(input))
    input
}


#' Change classes of multiple columns
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @param numeric_columns chr, vector of names of columns that will be 
#' converted to numeric.
#' @param character_columns chr, vector of names of colums taht will be 
#' converted to character.
#' @param factor_columns chr, vector of names of columns that will be 
#' converted to factor.
#' @return data.table
#' @keywords internal
.fixColumnTypes = function(input, numeric_columns = NULL, 
                           character_columns = NULL,
                           factor_columns = NULL) {
    # TODO: switch [[<- to [, .()]<-
    for (column in factor_columns) {
        input[[column]] = factor(input[[column]])
    }
    for (column in numeric_columns) {
        input[[column]] = as.numeric(as.character(input[[column]]))
    }
    for (column in character_columns) {
        input[[column]] = as.character(input[[column]])
    }
    input
}


#' Set column to a single value
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @param fill_list named list, names correspond to column names, elements 
#' to values that will be used in the columns.
#' @return data.table
#' @keywords internal
.fillValues = function(input, fill_list) {
    if (length(fill_list) > 0) {
        input[, names(fill_list) := as.list(unname(fill_list))]
    }
}


#' Change column names to match read.table/read.csv/read.delim conventions
#' @param col_names chr, vector of column names
#' @return character vector
#' @keywords internal
.standardizeColnames = function(col_names) {
    col_names = gsub(" ", ".", col_names, fixed = TRUE)
    col_names = gsub("\\[|\\]|\\%", ".", col_names, fixed = FALSE)
    col_names = gsub("/", "", col_names, fixed = TRUE)
    col_names = gsub("#", "X.", col_names, fixed = TRUE)
    gsub("[\\.]+", "", col_names)
}


#' Get intensity columns from wide-format data
#' @param col_names names of columns, where some of the columns store intensity
#' value for different channels
#' @param ... varying number of strings that define channel columns.
#' @return character vector of column names that correspond to channel intensities
#' @keywords internal
.getChannelColumns = function(col_names, ...) {
    all_patterns = unlist(list(...))
    channel_filter = rep(TRUE, length(col_names))
    for (pattern in all_patterns) {
        channel_filter = channel_filter & grepl(pattern, col_names, fixed = TRUE)
    }
    col_names[channel_filter]
}
