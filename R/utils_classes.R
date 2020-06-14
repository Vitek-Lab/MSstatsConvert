setOldClass("data.frame")

setClass("MSstatsValidated", contains = "data.frame")
setOldClass("MSstatsValidated", S4Class = "MSstatsValidated")

#' @importFrom methods new
#' @keywords internal
.MSstatsFormat = function(input) {
    new("MSstatsValidated", as.data.frame(input))
}
