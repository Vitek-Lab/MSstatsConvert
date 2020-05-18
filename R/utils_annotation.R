#' Create annotation
#' @param input data.table preprocessed by one of the .cleanRaw* functions
#' @param annotation data.table 
#' @param ... key-value pairs, where keys are names of columns of `annotation` 
.makeAnnotation = function(input, annotation, ...) {
    all_columns = unlist(list(...))
    annotation = .getDataTable(annotation)
    colnames(annotation) = .standardizeColnames(colnames(annotation))
    if (length(all_columns) > 0) {
        annotation = .updateColnames(annotation, 
                                     unname(all_columns),
                                     names(all_columns))
    }
    joint_columns = intersect(colnames(input), colnames(annotation)) 
    if (all(colnames(annotation) %in% joint_columns)) {
        NULL
    } else {
        annotation
    }
    # TODO: checks
}


#' Merge annotation with feature data
#' @param data.table preprocessed by one of the .cleanRaw* functions.
#' @param annotation data.table with annotation
#' @return data.table
#' @keywords internal 
.mergeAnnotation = function(input, annotation) {
    if (!is.null(annotation)) {
        cols = intersect(colnames(input), colnames(annotation))
        input = merge(input, annotation, by = cols, all.x = TRUE)
    }
    if (any(is.na(input$Condition))) {
        msg = "Condition in the input file must match condition in annotation"
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    if (is.element("Channel", colnames(input))) {
        if (!all(unique(annotation$Channel) %in% unique(input$Channel)))
            msg = "Please check the annotation file. The channel name must be matched with that in input data "
            getOption("MSstatsLog")("ERROR", msg)
            stop(msg)
    }
    input
}