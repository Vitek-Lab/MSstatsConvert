setClass("MSstatsValidated", contains = "data.frame")

setClass("MSstatsLabelFree", contains = "MSstatsValidated")

setClass("MSstatsLabeled", contains = "MSstatsValidated")

setClass("MSstatsTMT", contains = "MSstatsValidated")

.MSstatsFormat = function(input) {
    if (is.element("Channel", colnames(input))) {
        new("MSstatsTMT", new("MSstatsValidated", as.data.frame(input)))
    } else {
        if (length(unique(input$IsotopeLabelType)) == 1L) {
            new("MSstatsLabelFree", new("MSstatsValidated", 
                                        as.data.frame(input)))
        } else {
            new("MSstatsLabeled", new("MSstatsValidated", 
                                      as.data.frame(input)))
        }
    }
}


is(.MSstatsFormat(input))
is.data.frame(.MSstatsFormat(input))
inherits(.MSstatsFormat(input), "data.frame")
