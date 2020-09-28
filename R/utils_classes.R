setOldClass("data.frame")

setClass("MSstatsValidated", contains = "data.frame")
setOldClass("MSstatsValidated", S4Class = "MSstatsValidated")

#' Output format for further analysis by MSstats
#' @importFrom methods new
#' @return object of class MSstatsValidated that inherits from data.frame
#' @keywords internal
.MSstatsFormat = function(input) {
    new("MSstatsValidated", as.data.frame(input))
}
