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
                msg = paste("** Please check the annotation file.",
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
        msg = "** Run annotation merged with quantification data."
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    if (any(is.na(input$Condition))) {
        msg = "** Condition in the input file must match condition in annotation."
        getOption("MSstatsLog")("WARN", msg)
        getOption("MSstatsMsg")("WARN", msg)
    }
    input
}


#' Check if the annotation is valid
#' @param annotation annotation created by the MSstatsMakeAnnotation function
#' @return TRUE invisibly if the annotation is correct, throws an error otherwise
#' @keywords internal
.checkAnnotation = function(annotation) {
    if (is.element("Channel", colnames(annotation))) {
        missing_cols = setdiff(
            c("Run", "TechRepMixture", "Fraction", "Mixture", 
              "Channel", "Condition", "BioReplicate"),
            colnames(annotation)
        )
    } else {
        missing_cols = setdiff(c("Run", "Condition", "BioReplicate"),
                               colnames(annotation))
        counts = annotation[, list(n_rows = .N), by = "Run"]
        if (any(counts$n_rows > 1)) {
            msg = paste("** Please check annotation.",
                        "Each MS run (Raw.file) can\'t have multiple",
                        "conditions or BioReplicates.")
            getOption("MSstatsLog")("ERROR", msg)
            stop(msg)
        }
    }
    if (length(missing_cols) > 0) {
        msg = paste("** Columns", paste(missing_cols, sep = ", ", 
                                        collapse = ", "),
                    "missing in the annotation.", 
                    "Please check the annotation file.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    invisible(TRUE)
}