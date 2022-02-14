setOldClass("data.frame")

setClass("MSstatsValidated", contains = "data.frame")
setOldClass("MSstatsValidated", S4Class = "MSstatsValidated")

#' Output format for further analysis by MSstats
#' @param input data.table
#' @importFrom methods new
#' @return object of class MSstatsValidated that inherits from data.frame
#' @keywords internal
.MSstatsFormat = function(input) {
    input = .selectMSstatsColumns(input)
    new("MSstatsValidated", as.data.frame(input))
}


#' Convert output of converters to data.frame
#' @param x object of class MSstatsValidated
#' @return data.frame
#' @export
as.data.frame.MSstatsValidated = function(x) {
  as.data.frame(unclass(x))
}

#' Convert output of converters to data.table
#' @param x object of class MSstatsValidated
#' @return data.tables
#' @export
as.data.table.MSstatsValidated = function(x) {
  data.table::as.data.table(unclass(x))
}