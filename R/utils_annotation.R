#' Create annotation
#' @param input data.table preprocessed by one of the .cleanRaw* functions
#' @param annotation data.table 
#' @param ... key-value pairs, where keys are names of columns of `annotation` 
#' @keywords internal
.makeAnnotation = function(input, annotation, ...) {
    if (!is.null(annotation)) {
        all_columns = unlist(list(...))
        annotation = .getDataTable(annotation)
        if (length(all_columns) > 0) {
            colnames(annotation) = .updateColnames(annotation, 
                                                   unname(all_columns),
                                                   names(all_columns))
        }
        if (is.element("Channel", colnames(annotation))) {
            annotation$Channel = .standardizeColnames(annotation$Channel)
        }
        annotation[, !duplicated(colnames(annotation)), with = FALSE]
    } else {
        NULL
    }
}


#' Merge annotation with feature data
#' @param data.table preprocessed by one of the .cleanRaw* functions.
#' @param annotation data.table with annotation
#' @return data.table
#' @keywords internal 
.mergeAnnotation = function(input, annotation) {
    if (!is.null(annotation)) {
        if (is.element("Channel", colnames(input))) {
            if (!all(unique(annotation$Channel) %in% unique(input$Channel))) {
                msg = "Please check the annotation file. The channel name must be matched with that in input data "
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        cols = intersect(colnames(input), colnames(annotation))
        cols = setdiff(cols, c("Run", "Channel"))
        annotation_cols = intersect(c("Run", "Channel"), colnames(input))
        input = merge(input[, !(colnames(input) %in% cols), with = FALSE], 
                      annotation, by = annotation_cols, all.x = TRUE, sort = FALSE)
    }
    if (any(is.na(input$Condition))) {
        msg = "Condition in the input file must match condition in annotation"
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    input
}
