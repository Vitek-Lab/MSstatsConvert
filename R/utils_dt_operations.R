#' Read file from a provided path or convert given data.frame to data.table
#' 
#' @importFrom data.table as.data.table fread
#' @keywords internal
.getDataTable = function(input, ...) {
    if (is.data.frame(input)) {
        as.data.table(input)
    } else {
        data.table::fread(input, showProgress = FALSE, ...)
    }
}


#' Update specified column names while retaining the rest
#' @param input data.table preprocessed by one of the `cleanRaw*` functions.
#' @param old_names chr, vector of column names that will be changed.
#' @param new_names chr, vector of new columns names in order corresponding 
#' to that in `old_names` parameter.
#' @return character vector
#' @keywords internal
.updateColnames = function(input, old_names, new_names) {
    column_update = new_names
    names(column_update) = old_names
    columns <- colnames(input)
    not_changing <- setdiff(columns, names(column_update))
    column_update[not_changing] <- not_changing
    unname(column_update[columns])
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
#' @param fill_vector named vector, names correspond to column names, elements 
#' to values that will be used in the columns.
#' @return data.table
#' @keywords internal
.fillValues = function(input, fill_vector) {
    for (column in names(fill_vector)) {
        input[[column]] = fill_vector[column]
    }
    input
}


#' Change column names to match read.table/read.csv/read.delim conventions
#' @param col_names chr, vector of column names
#' @return character vector
#' @keywords internal
.standardizeColnames = function(col_names) {
    col_names = gsub(" ", ".", col_names, fixed = TRUE)
    gsub("#", "X.", col_names, fixed = TRUE)
}
