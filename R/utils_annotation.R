.getChannelColumns = function(col_names, ...) {
    all_patterns = unlist(list(...))
    channel_filter = rep(TRUE, length(col_names))
    for (pattern in all_patterns) {
        channel_filter = channel_filter & grepl(pattern, col_names, fixed = TRUE)
    }
    col_names[channel_filter]
}

# COMMON DATA PREPROCESSING
.checkColumnsGen = function(df_name, expected_columns, actual_columns, 
                            lgl_fun, msg_sep) {
    if (!lgl_fun(expected_columns %in% actual_columns)) {
        missing_columns = setdiff(expected_columns, actual_columns)
        stop(paste("Missing columns in", paste0(df_name, ":"), 
                   paste(missing_columns, sep = msg_sep, collapse = msg_sep)))
    }
}

.checkColumns = function(df_name, expected_columns, 
                         actual_columns, type = "required") {
    if (type == "required") {
        .checkColumnsGen(df_name, expected_columns, actual_columns, all, ", ")    
    } else {
        .checkColumnsGen(df_name, expected_columns, actual_columns, any, " or ")        
    }
    
}


.pickAnnotation = function(annotation, backup_annotation, columns_definition,
                           backup_columns_definition) {
    # Account for the case when user provides a path to file.
    if (is.null(annotation)) {
        .checkColumns("Annotation", names(backup_columns_definition), 
                      colnames(backup_annotation))
        list(df = backup_annotation,
             cols = backup_columns_definition)
    } else {
        .checkColumns("Annotation", names(columns_definition), 
                      colnames(annotation))
        list(df = annotation,
             cols = columns_definition)
    }
}

.checkAnnotationValidity = function(annotation) {
    counts_in_run = xtabs(~ Run, as.data.frame(annotation))
    if (any(counts_in_run > 1)) {
        stop('Please check annotation. Each MS Run must have a single condition and biological replicate')
    }
}

#' @importFrom data.table fread as.data.table
.makeAnnotation = function(annotation_source, columns_definition,
                           backup_annotation_source = NULL,
                           backup_columns_definition = NULL) {
    if (is.null(annotation_source) & is.null(backup_annotation_source)) {
        stop("Please provide annotation information")
    }
    if (is.null(backup_columns_definition) & !is.null(backup_annotation_source)) {
        backup_columns_definition = columns_definition
    }
    if (is.character(annotation_source)) {
        annotation_source = fread(annotation_source)
    }
    if (is.character(backup_annotation_source)) {
        backup_annotation_source = fread(backup_annotation_source)
    }
    
    annotation_list = .pickAnnotation(annotation_source, backup_annotation_source,
                                      columns_definition, backup_columns_definition)
    annotation_list[["df"]] = as.data.table(annotation_list[["df"]])
    colnames(annotation_list[["df"]]) = .updateColnames(annotation_list[["df"]], 
                                                        names(annotation_list[["cols"]]),
                                                        annotation_list[["cols"]])
    annotation_list[["df"]] = unique(annotation_list[["df"]][, annotation_list[["cols"]], with = FALSE])
    .checkAnnotationValidity(annotation_list[["df"]])
    annotation_list[["df"]]
}

.mergeAnnotation = function(input, annotation) {
    cols = intersect(colnames(input), colnames(annotation))
    input = merge(input, annotation, by = cols, all.x = TRUE)
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