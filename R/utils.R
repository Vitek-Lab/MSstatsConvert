#' A dummy function to store shared documentation items.
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

.updateColnames = function(data_frame, column_update) {
    columns <- colnames(data_frame)
    not_changing <- setdiff(columns, names(column_update))
    column_update[not_changing] <- not_changing
    unname(column_update[columns])
}


# ALL DATA CHECKS!!
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
    invisible(TRUE)
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
    data_frame[, !(colnames(data_frame) %in% columns_to_remove)]
}

.fixColumnTypes = function(data_frame, numeric_columns, character_columns,
                           factor_columns) {
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
    if(is.null(annotation)) {
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

.checkAnnotationValidity = function(annotation){
    counts_in_run = xtabs(~ Run, annotation)
    if (any(counts_in_run > 1)) {
        stop('Please check annotation. Each MS Run must have a single condition and biological replicate')
    }
}

.makeAnnotation = function(annotation_source, columns_definition,
                           backup_annotation_source = NULL,
                           backup_columns_definition = NULL) {
    if(is.null(annotation_source) & is.null(backup_annotation_source)) {
        stop("Please provide annotation information")
    }
    if(is.null(backup_columns_definition) & !is.null(backup_annotation_source)) {
        backup_columns_definition = columns_definition
    }
    annotation_list = .pickAnnotation(annotation_source, backup_annotation_source,
                                      columns_definition, backup_columns_definition)
    colnames(annotation_list[["df"]]) = .updateColnames(annotation_list[["df"]], 
                                                        annotation_list[["cols"]])
    annotation_list[["df"]] = unique(annotation_list[["df"]][, annotation_list[["cols"]]])
    .checkAnnotationValidity(annotation_list[["df"]])
    annotation_list[["df"]]
}

.findAvailable = function(possibilities, option_set) {
    chosen = option_set[option_set %in% possibilities]
    if(length(chosen) == 0) {
        NULL 
    } else {
        chosen
    }
}

.removeSharedPeptides = function(data_frame, proteins_column, peptides_column) {
    unique_pairs = unique(data_frame[, c(proteins_column, peptides_column)])
    protein_counts = aggregate(x = unique_pairs[[proteins_column]], 
                               by = list(peptide = unique_pairs[[peptides_column]]),
                               length)
    counts = protein_counts[["x"]]
    names(counts) = protein_counts[["peptide"]]
    if(length(counts) == 0) {
        data_frame
    } else {
        data_frame[counts[data_frame[[peptides_column]]] == 1, ]    
    }
    # TODO: message for the user / log
}

.handleSharedPeptides = function(data_frame, proteins_column, peptides_column,
                                 remove_shared = TRUE) {
    if(remove_shared) {
        .removeSharedPeptides(data_frame, proteins_column, peptides_column)
    } else {
        data_frame
    }
}

.filterByScore = function(data_frame, score_column, score_threshold, direction,
                          behavior, fill_value = NULL) {
    if(direction == "greater") {
        score_filter = data_frame[[score_column]] > score_threshold
    } else {
        score_filter = data_frame[[score_column]] < score_threshold
    }
    score_filter = score_filter & !is.na(data_frame[[score_column]])
    if(behavior == "remove") {
        data_frame[score_filter, ]    
    } else {
        data_frame[score_filter, ] = fill_value
        data_frame
    }
}

.handleFiltering = function(data_frame, score_column, score_threshold, 
                            direction, behavior, fill_value = NULL, 
                            drop_column = TRUE, filter = TRUE) {
    if(filter) {
        result = .filterByScore(data_frame, score_column, score_threshold, behavior,
                                fill_value)
    } else {
        result = data_frame
    }
    if(drop_column) {
        .removeColumns(result, score_column)
    } else {
        result
    }
}

.fillValues = function(data_frame, fill_vector) {
    for(column in names(fill_vector)) {
        data_frame[[column]] = fill_vector[column]
    }
    data_frame
}


.makeFeature = function(data_frame, feature_columns) {
    
}

.filterMissingFeatures = function(data_frame, feature) {
    
}

.filterFewMeasurements = function(data_frame, features) {
    
}

.summarizeMultipleMeasurements = function(data_frame, features) {
    
}

.cleanByFeature = function(data_frame, feature_columns, summarize_function,
                           handle_few_measurements,                         
) {
 # 1. Make features
    features = .makeFeature(data_frame, feature_columns)
 # 2. Remove features that are all 0/NA
    filtered = .filterMissingFeatures(data_frame, features)
 # 3. Remove features with few measurements
    filtered = .filterFewMeasurements(data_frame, features)
 # 4. Summarize multiple measurements per feature
    .summarizeMultipleMeasurements(data_frame, features)
}

# DIA-Umpire
# # 4 files: 
# # - fragment ion summary
# # - peptide summary
# # - protein summary
# MaxQuant
# # 2 or 3 files:
# # - evidence.txt (feature level data)
# # - annotation.txt (Raw.file, Condition, BioReplicate, Run, IsotopeLabelType)
# # - proteinGroups.txt (protein group ID). We can use "Proteins" column from evidence instead
# OpenMS
# # 1 or 2 files:
# # - feature-level data
# # - annotation (Condition, BioReplicate, Run). Can be replaced with columns from the feature-level data 
# OpenSWATH
# # 2 files:
# # - feature-level data
# # - annotation file ('Condition', 'BioReplicate', 'Run')
# Progenesis
# # - feature-level data ('Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required)
# # - annotation file (Condition, BioReplicate, Run)
# Skyline
# # 1 or 2 files:
# # - feature-level data
# # - annotation file - can be replaced with columns "Run", "Condition", 'BioReplicate' in the feature-level data
# Spectronaut
# # 1 or 2 files:
# # - feature-level data
# # - annotation file - can be replaced with "R.FileName", "R.Condition", "R.Replicate" in the feature-level data
# Proteome Discoverer
# # 2 files:
# # - feature-level data 
# # - annotation file - 'Condition', 'BioReplicate', 'Run'
