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

.checkAnnotationValidity = function(annotation) {
    counts_in_run = xtabs(~ Run, as.data.frame(annotation))
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
        .removeSharedPeptides(data_frame, protein_column, peptide_column)
    } else {
        data_frame
    }
}

.handleOxidationPeptides = function(data_frame, sequence_column, 
                                    oxidation_pattern, remove) {
    if(remove) {
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

.handleFiltering = function(data_frame, score_column, score_threshold, 
                            direction, behavior, fill_value = NULL, 
                            drop_column = TRUE, filter = TRUE) {
    if(filter) { 
        result = .filterByScore(data_frame, score_column, score_threshold, 
                                direction, behavior, fill_value)
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

.filterSmallIntensities = function(data_frame, threshold) {
    threshold_filter = data_frame[["Intensity"]] > threshold
    threshold_filter = threshold_filter & !is.na(data_frame[["Intensity"]])
    threshold_filter
}

.makeFeatures = function(data_frame, feature_columns) {
    gsub(" " , "", apply(data_frame[, feature_columns, drop = FALSE,
                                    with = FALSE], 1, 
                         function(x) paste(x, sep = "_", collapse = "_")))
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
    if(remove_single_feature & length(single_feature) > 0) {
        data_frame = data_frame[!(data_frame[["ProteinName"]] %in% single_feature), ]
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
    data_frame = .summarizeMultipleMeasurements(data_frame, summarize_function)
    data_frame = .filterFewMeasurements(data_frame, 0, handle_few_measurements)
    data_frame
}

.handleDecoyProteins = function(data_frame, decoy_column, decoy_symbols, drop = TRUE) {
    decoy_index = which(colnames(data_frame) == decoy_column)
    if(drop) {
        data_frame[!(data_frame[[decoy_column]] %in% decoy_symbols), -decoy_index]    
    } else {
        data_frame[!(data_frame[[decoy_column]] %in% decoy_symbols), ]
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

.cleanRawOpenSWATH = function(openswath_input) {
    os_cols = c("ProteinName", "FullPeptideName", "Charge", "filename", 
                "aggr_Fragment_Annotation", "aggr_Peak_Area")
    input = data.table::as.data.table(openswath_input[, os_cols])
    colnames(input) = .updateColnames(
        input,
        c("FullPeptideName", "Charge", "filename"),
        c("PeptideSequence", "PrecursorCharge", "Run"))
    input = input[, lapply(.(aggr_Fragment_Annotation, aggr_Peak_Area), 
                           function(x) unlist(tstrsplit(x, ";", fixed = TRUE))),
                  by = .(ProteinName, PeptideSequence, PrecursorCharge, Run)]
    colnames(input) = .updateColnames(input, c("V1", "V2"), 
                                      c("FragmentIon", "Intensity"))
    input[, c("PeptideSequence", "FragmentIon")] = input[, lapply(.(PeptideSequence, FragmentIon), 
                                                                  function(x) gsub(":", "_", x))]
    input[["Intensity"]] = as.numeric(input[["Intensity"]])
    input[input[["Intensity"]] < 1, "Intensity"] = NA
    input = .fillValues(input, c("ProductCharge" = NA, "IsotopeLabelType" = "L"))
    input
}

.cleanRawPD = function(pd_input, quantification_column, proteinID_column,
                       sequence_column, filter_num_col) {
    colnames(pd_input) = .standardizeColnames(pd_input)
    which.quantification = .findAvailable(c("Intensity", "Area"),
                                          colnames(pd_input),
                                          "Precursor.Area")
    which.quantification = .isLegalValue(which.quantification, 
                                         legal_values = c("Intensity", "Area", "Precursor.Area",
                                                          "Precursor.Abundance"),
                                         message = "Please select a column to be used for quantified intensities among four options: ")
    
    which.proteinid = .findAvailable(c("Protein.Accessions", 
                                       "Master.Protein.Accessions"),
                                     colnames(pd_input),
                                     "Protein.Group.Accessions")
    which.proteinid = .isLegalValue(which.proteinid, 
                                    legal_values = c("ProteinAccessions", 
                                                     "Master.Protein.Accessions",
                                                     "Protein.Group.Accessions"),
                                    message = "Please select a column to be used as protein IDs among three options: ")
    which.sequence = .findAvailable("Annotated.Sequence", colnames(pd_input), 
                                    "Sequence")
    which.sequence = .isLegalValue(which.sequence, 
                                   legal_values = c("Annotated.Sequence", "Sequence"),
                                   message = "Please select peptide sequence column between two options: ")
    
    if(filter_num_col) {
        pd_input = pd_input[pd_input[["X..Proteins"]] == '1', ]
    }
    pd_cols = c(proteinID_column, "X..Proteins", sequence_column, 
                "Modifications", "Charge", "Spectrum.File", quantification_column)
    if (any(is.element(colnames(pd_input), "Fraction"))) {
        pd_cols = c(pd_cols, "Fraction")
    }
    input = data.table::as.data.table(pd_input[, pd_cols])
    colnames(input) = .updateColnames(
        input,
        c(proteinID_column, sequence_column, "Spectrum.File", quantification_column),
        c("ProteinName", "PeptideSequence", "Run", "Intensity"))
    input[["PeptideModifiedSequence"]] = paste(input[["PeptideSequence"]], 
                                               input[["Modifications"]], 
                                               sep = "_")
    input[, -which(colnames(input) == "PeptideSequence")]
}

.standardizeColnames = function(col_names) {
    gsub(" ", ".", col_names, fixed = TRUE)
}


.cleanRawProgenesis = function(prog_input, runs, fix_colnames = TRUE) {
    prog_input = data.table::as.data.table(prog_input)
    colnames(prog_input) = .standardizeColnames(colnames(prog_input))
    if(fix_colnames) {
        prog_input = prog_input[-1, ]
        colnames(prog_input) = prog_input[1, ]
        prog_input = prog_input[-1, ]
    }
    protein_col = .findAvailable(c("Protein", "Accession"), 
                                 colnames(prog_input))
    colnames(prog_input) = .updateColnames(prog_input, 
                                           c(protein_col, "Charge"), 
                                           c("ProteinName", "PrecursorCharge"))
    
    nonmissing_prot = !is.na(prog_input[["ProteinName"]]) & prog_input[["ProteinName"]] != ""
    nonmissing_pept = !is.na(prog_input[["Sequence"]]) & prog_input[["Sequence"]] != ""
    prog_input = prog_input[nonmissing_prot & nonmissing_pept, ]
    prog_input[["PeptideModifiedSequence"]] = paste(input[["Sequence"]],
                                                    input[["Modifications"]],
                                                    sep = "")
    prog_input <- prog_input[!duplicated(prog_input), ] # dubious performance-wise
    if(is.element("Use.in.quantitation", colnames(prog_input))) {
        prog_input = prog_input[prog_input[["Use.in.quantitation"]], ]
        # TODO: consider character version if these files ever import this 
        # columns as character
        prog_input = prog_input[, Use.in.quantitation := NULL]
    }
    feature_cols = c("ProteinName", "PeptideModifiedSequence",
                     "PrecursorCharge", "Fraction")
    prog_cols = intersect(colnames(prog_input), 
                          c(feature_cols, runs))
    prog_input = prog_input[, prog_cols]
    prog_input = data.table::melt(prog_input,
                                  id.vars = feature_cols,
                                  measure_vars = runs,
                                  variable.name = "Run",
                                  value.name = "Intensity")
    prog_input[["Intensity"]] = as.numeric() # Remove if not needed
    prog_input
}

.cleanRawSpectronaut = function(spec_input) {
    spec_input = data.table::as.data.table(spec_input)
    colnames(spec_input) = .standardizeColnames(colnames(spec_input))
    
    spec_input = spec_input[spec_input[["F.FrgLossType"]] == "noloss", ]
    spec_input = spec_input[!spec_input[["F.ExcludedFromQuantification"]], ]
    # XIC quality. TODO: explain in documentation
    
    f_charge_col = .findAvailable(c("F.Charge", "F.FrgZ"), colnames(spec_input))
    pg_qval_col = .findAvailable(c("PG.Qvalue"), colnames(spec_input))
    spec_input = .selectColumns(
        spec_input, 
        c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon", 
          f_charge_col, "R.FileName", "EG.Qvalue", pg_qval_col, paste0("F.", intensity)))
    colnames(spec_input) = .updateColnames(
        spec_input, 
        c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon",
          f_charge_col, "R.FileName", "EG.Qvalue", paste0("F.", intensity)),
        c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
          "ProductCharge", "Run", "Qvalue", "Intensity"))
    spec_input
}

.aggregateMonoisotopicPeaks <- function(data_frame) {
    data_frame[["pepprecursor"]] <- paste(data_frame[["PeptideSequence"]], 
                                          data_frame[["PrecursorCharge"]], 
                                          sep = "_")
    data_frame <- data_frame[!is.na(data_frame[["Intensity"]]), ]
    standard.info <- unique(data_frame[, c("ProteinName", "PeptideSequence", 
                                           "PrecursorCharge", "StandardType")])
    data_w <- data.table::dcast(data = data_frame, pepprecursor ~ Run, 
                                value.var = "Intensity", 
                                fun.aggregate = function(x) sum(x, na.rm = TRUE), 
                                fill = NA_real_) 
    newdata <- data.table::melt(data_w, id.vars = c("pepprecursor"))
    colnames(newdata) = .updateColnames(newdata, c("variable", "value"),
                                        c("Run", "Intensity"))
    uniinfo <- unique(data_frame[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "pepprecursor")])	
    data_frame <- merge(newdata, uniinfo, by = "pepprecursor")
    data_frame <- merge(data_frame, standard.info, 
                        by = c("ProteinName", "PeptideSequence", "PrecursorCharge"))
    data_frame = .fillValues(data_frame, c("FragmentIon" = "sum",
                                           "ProductCharge" = NA,
                                           "IsotopeLabelType" = "L"))
    data_frame = data_frame[, pepprecursor := NULL]
    data_frame
}

.cleanRawSkyline = function(sl_input) {
    colnames(sl_input) = gsub("\\.", "", colnames(sl_input))
    colnames(sl_input) = .updateColnames(sl_input, c("FileName", "Area"),
                                         c("Run", "Intensity"))
    colnames(sl_input) = .standardizeColnames(sl_input)
    
    sl_input = data.table::as.data.table(sl_input) 
    if(is.element("PeptideSequence", colnames(sl_input))) {
        sl_input = sl_input[, .(PeptideSequence) := NULL]
    }
    colnames(sl_input) = .updateColnames(sl_input, "PeptideModifiedSequence",
                                         "PeptideSequence")
    sl_input[["Intensity"]] = as.numeric(sl_input[["Intensity"]])
    if(is.element("DetectionQValue", colnames(sl_input))) {
        sl_input[["DetectionQValue"]] = as.numeric(as.character(sl_input[["DetectionQValue"]]))    
    }
    if(is.character(sl_input[["Truncated"]])) {
        sl_input[["Truncated"]] = sl_input[["Truncated"]] == "True"
    }
    sl_input[["Truncated"]] = as.integer(sl_input[["Truncated"]])
    
    sl_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                "FragmentIon", "ProductCharge", "IsotopeLabelType", "Condition",
                "BioReplicate", "Run", "Intensity", "StandardType")
    sl_cols = c(sl_cols, "Fraction", "DetectionQValue", "Truncated")
    sl_input = sl_input[, intersect(sl_cols, colnames(sl_input)), with = FALSE]
    sl_input
}

.cleanRawOpenMS = function(om_input) {
    colnames(om_input) = .standardizeColnames(om_input)
    om_input = data.table::as.data.table(om_input)
    om_input[["Intensity"]] = as.numeric(om_input[["Intensity"]])
    if(!is.element("IsotopeLabelType", colnames(om_input))) {
        om_input = .fillValues(om_input, c("IsotopeLabelType" = "L"))
    }
    om_input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                 "FragmentIon", "ProductCharge", "IsotopeLabelType",
                 "Condition", "BioReplicate", "Run", "Intensity"),
             with = FALSE]
}
