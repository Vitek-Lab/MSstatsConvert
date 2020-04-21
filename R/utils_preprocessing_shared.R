#' Read file from a provided path or convert given data.frame to data.table
#' 
#' @importFrom data.table as.data.table fread
#' @keywords internal
.getDataTable = function(input, ...) {
    if (is.data.frame(input)) {
        as.data.table(input)
    } else {
        data.table::fread(input, showProgress = FALSE, ...)
    }
}

#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' 
#' @param fewMeasurements 'remove'(default) will remove the features that have 1 or 2 measurements across runs.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.
#' @param removeMpeptides TRUE will remove the peptides including 'M' sequence. FALSE is default.
#' 
#' @keywords internal
#' 

.documentFunction = function(fewMeasurements, 
                             useUniquePeptide,
                             summaryforMultipleRows, 
                             removeProtein_with1Feature,
                             removeProtein_with1Protein,
                             removeOxidationMpeptides,
                             removeMpeptides) {
    
}

.updateColnames = function(data_frame, old_names, new_names) {
    column_update = new_names
    names(column_update) = old_names
    columns <- colnames(data_frame)
    not_changing <- setdiff(columns, names(column_update))
    column_update[not_changing] <- not_changing
    unname(column_update[columns])
}


# TODO: ALL DATA CHECKS!!
.makeInclusionErrorMessage = function(parameter_name, legal_values, 
                                      information) {
    if(is.null(information)) {
        paste("Parameter", parameter_name, "must be one of", 
              paste(legal_values, sep = ", ", collapse = ", "))
    } else {
        paste(information, paste(legal_values, sep = ", ", collapse = ", "))
    }
}

.isLegalValue = function(parameter, legal_values = NULL, 
                         can_be_null = FALSE, information = NULL) {
    parameter_name = deparse(substitute(parameter))
    if(is.null(parameter)) {
        if(!can_be_null) {
            stop(paste("Parameter", parameter_name, "cannot be NULL"))    
        }
    } else {
        if(!is.null(legal_values)) {
            if(!is.element(parameter, legal_values)) {
                stop(.makeInclusionErrorMessage(parameter_name, legal_values,
                                                information))
            }    
        }
    }
    parameter
}

# COMMON DATA PREPROCESSING
.checkColumnsGen = function(df_name, expected_columns, actual_columns, 
                            lgl_fun, msg_sep) {
    if(!lgl_fun(expected_columns %in% actual_columns)) {
        missing_columns = setdiff(expected_columns, actual_columns)
        stop(paste("Missing columns in", paste0(df_name, ":"), 
                   paste(missing_columns, sep = msg_sep, collapse = msg_sep)))
    }
}

.checkColumns = function(df_name, expected_columns, 
                         actual_columns, type = "required") {
    if(type == "required") {
        .checkColumnsGen(df_name, expected_columns, actual_columns, all, ", ")    
    } else {
        .checkColumnsGen(df_name, expected_columns, actual_columns, any, " or ")        
    }
    
}

.selectColumns = function(data_frame, column_names, label = "Input") {
    .checkColumns(label, column_names, colnames(data_frame))
    data_frame[, column_names]
}

.removeColumns = function(data_frame, columns_to_remove) {
    data_frame[, !(colnames(data_frame) %in% columns_to_remove), with = FALSE]
}

.fixColumnTypes = function(data_frame, numeric_columns = NULL, 
                           character_columns = NULL,
                           factor_columns = NULL) {
    for(column in factor_columns) {
        data_frame[[column]] = factor(data_frame[[column]])
    }
    for(column in numeric_columns) {
        data_frame[[column]] = as.numeric(as.character(data_frame[[column]]))
    }
    for(column in character_columns) {
        data_frame[[column]] = as.character(data_frame[[column]])
    }
    data_frame
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
    if(is.null(annotation_source) & is.null(backup_annotation_source)) {
        stop("Please provide annotation information")
    }
    if(is.null(backup_columns_definition) & !is.null(backup_annotation_source)) {
        backup_columns_definition = columns_definition
    }
    if(is.character(annotation_source)) {
        annotation_source = fread(annotation_source)
    }
    if(is.character(backup_annotation_source)) {
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

.findAvailable = function(possibilities, option_set, fall_back = NULL) {
    chosen = option_set[option_set %in% possibilities]
    if(length(chosen) != 1L) {
        if(is.null(fall_back)) {
            NULL
        } else {
            if(fall_back %in% possibilities) {
                fall_back
            } else {
                NULL
            }
        }
    } else {
        chosen
    }
}

.removeSharedPeptides = function(data_frame, proteins_column, peptides_column) {
    unique_pairs = unique(data_frame[, c(proteins_column, peptides_column), with = FALSE])
    protein_counts = aggregate(x = unique_pairs[[proteins_column]], 
                               by = list(peptide = unique_pairs[[peptides_column]]),
                               length)
    counts = protein_counts[["x"]]
    names(counts) = protein_counts[["peptide"]]
    if(length(counts) == 0) {
        data_frame
    } else {
        data_frame[counts[data_frame[[peptides_column]]] == 1L, ]    
    }
    # TODO: message for the user / log
}

.handleSharedPeptides = function(data_frame, remove_shared = TRUE,
                                 protein_column = "ProteinName",
                                 peptide_column = "PeptideSequence") {
    if(remove_shared) {
        getOption("MSstatsLog")("INFO", "Shared peptides are removed")
        getOption("MSstatsMsg")("INFO", "Shared peptides are removed")
        .removeSharedPeptides(data_frame, protein_column, peptide_column)
    } else {
        data_frame
    }
}

.handleOxidationPeptides = function(data_frame, sequence_column, 
                                    oxidation_pattern, remove) {
    if(remove) {
        msg = paste("Peptides containing", oxidation_pattern, "are removed")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        data_frame[!(grepl(oxidation_pattern, data_frame[[sequence_column]])), ]   
    } else {
        data_frame
    }
}

.filterByScore = function(data_frame, score_column, score_threshold, direction,
                          behavior, fill_value = NULL) {
    if(direction == "greater") {
        score_filter = data_frame[[score_column]] >= score_threshold
    } else {
        score_filter = data_frame[[score_column]] <= score_threshold
    }
    score_filter = score_filter & !is.na(data_frame[[score_column]])
    if(behavior == "remove") {
        data_frame[score_filter, ]    
    } else {
        data_frame[!score_filter, c("Intensity")] = fill_value
        data_frame
    }
}

.handleFiltering = function(input, score_column, score_threshold, 
                            direction, behavior, fill_value = NULL, 
                            drop_column = TRUE, filter = TRUE) {
    if (filter) { 
        result = .filterByScore(input, score_column, score_threshold, 
                                direction, behavior, fill_value)
    } else {
        result = input
    }
    
    if (drop_column) {
        .removeColumns(result, score_column)
    } else {
        result
    }
}


.fillValues = function(data_frame, fill_vector) {
    for (column in names(fill_vector)) {
        data_frame[[column]] = fill_vector[column]
    }
    data_frame
}

.filterSmallIntensities = function(data_frame, threshold) {
    threshold_filter = data_frame[["Intensity"]] > threshold
    threshold_filter = threshold_filter & !is.na(data_frame[["Intensity"]])
    threshold_filter
}

.makeFeatures = function(data_frame, feature_columns) {
    gsub(" ", "", do.call(.combine, dt[, feats, with = FALSE]))
}

.filterFewMeasurements = function(data_frame, min_intensity, handle_few) {
    int_filter = data_frame[["Intensity"]] > min_intensity
    int_filter = int_filter & !is.na(data_frame[["Intensity"]])
    counts = data_frame[int_filter, .(n_obs = length(Intensity)), 
                        by = .(feature)]
    if(handle_few == "remove") {
        not_few = unique(counts[["feature"]][counts[["n_obs"]] > 2])
    } else {
        not_few = unique(counts[["feature"]][counts[["n_obs"]] > 0])
    }
    data_frame[data_frame[["feature"]] %in% not_few, ]
    # TODO: compare performance to join
    # TODO: improve design by making minimum number a variable?
}

.summarizeMultipleMeasurements = function(data_frame, aggregator) {
    counts = data_frame[, ("n_obs" = length("Intensity")), by = ("feature"),
                        with = FALSE]
    if(any(counts[["n_obs"]] > length(unique(data_frame[["Run"]])))) {
        merge(data_frame[, .(Intensity = aggregator(Intensity)), 
                         by = list(Run, feature)],
              data_frame[, Intensity := NULL],
              by = c("Run", "feature")
              
        )
    } else {
        data_frame
    }
}

.handleSingleFeaturePerProtein = function(data_frame, remove_single_feature) {
    counts = data_frame[, .(n_obs = length(Intensity)), by = .(feature)]
    single_feature = counts[["feature"]][counts[["n_obs"]] <= 1]
    if (remove_single_feature & length(single_feature) > 0) {
        data_frame = data_frame[!(data_frame[["ProteinName"]] %in% single_feature), ]
        getOption("MSstatsLog")("INFO", "Proteins with a single feature are removed")
        getOption("MSstatsMsg")("INFO", "Proteins with a single feature are removed")
    } else {
        data_frame = data_frame
    }
    # TODO: message + logs
    feature_col = colnames(data_frame) == "feature"
    data_frame[, !feature_col, with = FALSE]
}


.cleanByFeature = function(data_frame, feature_columns, summarize_function,
                           handle_few_measurements) {
    data_frame[["feature"]] = .makeFeatures(data_frame, feature_columns)    
    data_frame = .filterFewMeasurements(data_frame, 1, "keep")
    getOption("MSstatsLog")("INFO", "Features with 1 or two measurements across runs are removed")
    getOption("MSstatsMsg")("INFO", "Features with 1 or two measurements across runs are removed")
    data_frame = .summarizeMultipleMeasurements(data_frame, summarize_function)
    getOption("MSstatsLog")("INFO", "Multiple measurements per run are aggregated")
    getOption("MSstatsMsg")("INFO", "Multiple measurements per run are aggregated")
    data_frame = .filterFewMeasurements(data_frame, 0, handle_few_measurements)
    data_frame
}

.handleDecoyProteins = function(data_frame, decoy_column, decoy_symbols, 
                                drop = TRUE, filter = TRUE) {
    decoy_index = which(colnames(data_frame) == decoy_column)  
    if(filter) {
        decoy_filter = !(data_frame[[decoy_column]] %in% decoy_symbols)
    } else {
        decoy_filter = rep(TRUE, nrow(data_frame))
    }
    
    if(drop) {
        data_frame[decoy_filter, -decoy_index, with = FALSE]    
    } else {
        data_frame[decoy_filter, ]
    }
}

.checkDDA = function(input) {
    # For now, assume Skyline input. Might need to be more general in the future
    fragment_ions = as.character(unique(input[["FragmentIon"]]))
    check_DDA = setdiff(c("precursor", "precursor [M+1]", "precursor [M+2]"), 
                        fragment_ions)
    frags = setdiff(fragment_ions, 
                    c('precursor', 'precursor [M+1]', 'precursor [M+2]'))
    precursors = intersect(fragment_ions, 
                           c("precursor", "precursor [M+1]", "precursor [M+2]"))
    ## if there are fragment ion and also have any 'precursor', it is the issue.
    if (length(frags) > 0 & length(precursors) > 0) {
        stop("** Please check precursors information. If your experiment is DIA, please remove the precursors. If your experiments is DDA, please check the precursor information.")
    }
    length(check_DDA) < 3
}

.filterManyColumns = function(data_frame, filter_columns, filter_symbols) {
    for(col in filter_columns) {
        has_col = is.element(col, colnames(data_frame))
        data_frame = .handleDecoyProteins(data_frame, col, filter_symbols, 
                                          has_col, has_col)        
    }
    data_frame
}
