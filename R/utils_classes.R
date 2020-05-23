setOldClass("data.frame")

setClass("MSstatsValidated", contains = "data.frame")
setOldClass("MSstatsValidated", S4Class = "MSstatsValidated")

setClass("MSstatsLabelFree", contains = "MSstatsValidated")
setOldClass("MSstatsLabelFree", S4Class = "MSstatsLabelFree")

setClass("MSstatsLabeled", contains = "MSstatsValidated")
setClass("MSstatsTMT", contains = "MSstatsValidated")


#' @importFrom methods new
#' @keywords internal
.MSstatsFormat = function(input) {
    if (is.element("Channel", colnames(input))) {
        new("MSstatsTMT", new("MSstatsValidated", as.data.frame(input)))
    } else {
        if (is.element("IsotopeLabelType", colnames(input))) {
            if (length(unique(input$IsotopeLabelType)) == 1L) {
                new("MSstatsLabelFree", new("MSstatsValidated", 
                                            as.data.frame(input)))
            } else {
                new("MSstatsLabeled", new("MSstatsValidated", 
                                          as.data.frame(input)))
            }
        }
        else {
            msg = "Incorrect data format: missing isotope label information"
            getOption("MSstatsLog")("ERROR", msg)
            stop(msg)
        }
    }
}
