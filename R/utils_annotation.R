#' Merge annotation with feature data
#' @param data.table preprocessed by one of the .cleanRaw* functions.
#' @param annotation data.table with annotation
#' @return data.table
#' @keywords internal 
.mergeAnnotation = function(input, annotation) {
    input$Run = .standardizeColnames(input$Run)
    if (!is.null(annotation)) {
        if (is.element("Channel", colnames(input))) {
            input$Channel = .standardizeColnames(input$Channel)
            if (!all(unique(annotation$Channel) %in% unique(input$Channel))) {
                msg = paste("Please check the annotation file.",
                            "The channel name must be matched",
                            "with that in input data.")
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        cols = intersect(colnames(input), colnames(annotation))
        cols = setdiff(cols, c("Run", "Channel"))
        annotation_cols = intersect(c("Run", "Channel"), colnames(input))
        input = merge(input[, !(colnames(input) %in% cols), with = FALSE], 
                      annotation, 
                      by = annotation_cols, all.x = TRUE, sort = FALSE)
    }
    if (any(is.na(input$Condition))) {
        msg = "Condition in the input file must match condition in annotation."
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    input
}
