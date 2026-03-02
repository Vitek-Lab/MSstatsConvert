#' Check validity of parameters to the `MSstatsImport` function.
#' @inheritParams MSstatsImport
#' @return none, throws an error if any of the assertions fail
#' @keywords internal
.checkMSstatsParams = function(input, annotation, 
                               feature_columns,
                               remove_shared_peptides,
                               remove_single_feature_proteins,
                               feature_cleaning) {
    checkmate::assertDataTable(input)
    checkmate::assertTRUE(is.character(annotation) | 
                              inherits(annotation, "data.frame") | 
                              is.null(annotation))
    checkmate::assertCharacter(feature_columns, min.len = 1)
    checkmate::assertLogical(remove_shared_peptides)
    checkmate::assertLogical(
        feature_cleaning[["remove_features_with_few_measurements"]], 
        )
    checkmate::assertLogical(remove_single_feature_proteins)
    checkmate::assertLogical(feature_cleaning[["remove_psms_with_all_missing"]],
                             null.ok = TRUE)
    checkmate::assertFunction(feature_cleaning[["summarize_multiple_psms"]])
}


#' Select an available options from a set of possibilities
#' @param possibilities possible legal values of a variable
#' @param option_set set of values that includes one of the `possibilities`
#' @param fall_back if there is none of the `possibilities` in `option_set`,
#' or there are multiple hits, default to `fall_back`
#' @return same as option_set, usually character
#' @keywords internal
.findAvailable = function(possibilities, option_set, fall_back = NULL) {
    chosen = option_set[option_set %in% possibilities]
    if (length(chosen) != 1L) {
        if (is.null(fall_back)) {
            NULL
        } else {
            if (fall_back %in% possibilities) {
                fall_back
            } else {
                NULL
            }
        }
    } else {
        chosen
    }
}


#' Generic parameter validation for all MSstats converters using configuration object
#' @param config A list containing all converter parameters. See details for required structure.
#' @details
#' The config list should contain the input and optionally other parameters:
#' - input: input data (required)
#' - annotation: annotation data (optional)
#' - intensity: intensity type (optional)
#' - filter_with_Qvalue: Q-value filter setting (default: FALSE)
#' - qvalue_cutoff: Q-value cutoff (default: 0.01)
#' - useUniquePeptide: unique peptide setting (default: TRUE)
#' - removeFewMeasurements: remove few measurements setting (default: TRUE)
#' - removeProtein_with1Feature: remove single feature proteins setting (default: FALSE)
#' - summaryforMultipleRows: aggregation function (default: max)
#' - calculateAnomalyScores: anomaly detection setting (default: FALSE)
#' - anomalyModelFeatures: anomaly model features (default: c())
#' - anomalyModelFeatureTemporal: temporal features (default: c())
#' - removeMissingFeatures: missing feature threshold (default: 0.5)
#' - anomalyModelFeatureCount: feature count for anomaly model (default: 100)
#' - runOrder: run order data (default: NULL)
#' - n_trees: number of trees (default: 100)
#' - max_depth: max tree depth (default: "auto")
#' - numberOfCores: number of cores (default: 1)
#' - use_log_file: logging setting (default: TRUE)
#' - append: append setting (default: FALSE)
#' - verbose: verbose setting (default: TRUE)
#' - log_file_path: log file path (default: NULL)
#' - excludedFromQuantificationFilter: filter setting (default: NULL)
#' @return NULL (throws error if validation fails)
#' @keywords internal
.validateMSstatsConverterParameters = function(config) {
    
    # Ensure config is a list
    if (!is.list(config)) {
        stop("Config must be a list")
    }
    
    # Define all defaults in one place
    defaults = list(
        input = NULL,
        annotation = NULL,
        intensity = NULL,
        excludedFromQuantificationFilter = NULL,
        filter_with_Qvalue = FALSE,
        qvalue_cutoff = 0.01,
        useUniquePeptide = TRUE,
        removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE,
        summaryforMultipleRows = max,
        calculateAnomalyScores = FALSE,
        anomalyModelFeatures = c(),
        anomalyModelFeatureTemporal = c(),
        removeMissingFeatures = 0.5,
        anomalyModelFeatureCount = 100,
        runOrder = NULL,
        n_trees = 100,
        max_depth = "auto",
        numberOfCores = 1,
        use_log_file = TRUE,
        append = FALSE,
        verbose = TRUE,
        log_file_path = NULL
    )
    
    # Merge config with defaults (config values override defaults)
    config = modifyList(defaults, config)
    
    # Input data validation (fail immediately if data is invalid)
    if (is.null(config$input)) {
        stop("Input data cannot be NULL")
    }
    
    if (is.character(config$input)) {
        if (!file.exists(config$input)) {
            stop("Input file does not exist: ", config$input)
        }
        # Quick file size check for very large files
        file_size_mb = file.size(config$input) / (1024^2)
        if (file_size_mb > 1000) {  # Warn for files > 1GB
            warning("Large input file detected (", round(file_size_mb, 1), " MB). ",
                    "Consider validating parameters on a subset first.")
        }
    } else if (is.data.frame(config$input) || data.table::is.data.table(config$input)) {
        # Quick structural validation
        if (nrow(config$input) == 0) {
            stop("Input data is empty (0 rows)")
        }
        if (ncol(config$input) == 0) {
            stop("Input data has no columns")
        }
    } else {
        stop("Input must be a file path, data.frame, or data.table")
    }
    
    # Annotation validation
    if (!is.null(config$annotation)) {
        if (is.character(config$annotation)) {
            if (!file.exists(config$annotation)) {
                stop("Annotation file does not exist: ", config$annotation)
            }
        } else if (!is.data.frame(config$annotation) && !data.table::is.data.table(config$annotation)) {
            stop("Annotation must be NULL, a file path, data.frame, or data.table")
        }
    }
    
    # Q-value filtering parameters
    checkmate::assertLogical(config$filter_with_Qvalue, len = 1)
    checkmate::assertNumber(config$qvalue_cutoff, lower = 0, upper = 1)
    
    # Common processing parameters
    checkmate::assertLogical(config$useUniquePeptide, len = 1)
    checkmate::assertLogical(config$removeFewMeasurements, len = 1)
    checkmate::assertLogical(config$removeProtein_with1Feature, len = 1)
    
    # Aggregation function validation
    if (!is.null(config$summaryforMultipleRows)) {
        checkmate::assertFunction(config$summaryforMultipleRows)
    }
    
    # Core system parameters
    checkmate::assertInt(config$numberOfCores, lower = 1)
    checkmate::assertLogical(config$use_log_file, len = 1)
    checkmate::assertLogical(config$append, len = 1)
    checkmate::assertLogical(config$verbose, len = 1)
    
    # Converter-specific boolean parameters (if provided)
    if (!is.null(config$excludedFromQuantificationFilter)) {
        checkmate::assertLogical(config$excludedFromQuantificationFilter, len = 1)
    }
    
    checkmate::assertLogical(config$calculateAnomalyScores, len = 1)
    
    # Anomaly detection parameter validation (converter-specific)
    if (config$calculateAnomalyScores) {
        # These validations only matter if anomaly detection is enabled
        if (length(config$anomalyModelFeatures) == 0) {
            stop("anomalyModelFeatures cannot be empty when calculateAnomalyScores=TRUE")
        }
        checkmate::assertCharacter(config$anomalyModelFeatures, min.len = 1)
        
        if (length(config$anomalyModelFeatureTemporal) > 0) {
            if (length(config$anomalyModelFeatureTemporal) != length(config$anomalyModelFeatures)) {
                stop("anomalyModelFeatureTemporal must have same length as anomalyModelFeatures or be empty")
            }
            valid_temporal = c("mean_decrease", "mean_increase", "dispersion_increase")
            invalid_temporal = config$anomalyModelFeatureTemporal[
                !is.null(config$anomalyModelFeatureTemporal) & 
                    !config$anomalyModelFeatureTemporal %in% valid_temporal
            ]
            if (length(invalid_temporal) > 0) {
                stop("Invalid temporal directions: ", paste(invalid_temporal, collapse = ", "),
                     ". Must be one of: ", paste(valid_temporal, collapse = ", "), " or NULL")
            }
        }
        
        checkmate::assertInt(config$n_trees, lower = 1)
        if (is.character(config$max_depth)) {
            checkmate::assertChoice(config$max_depth, choices = "auto")
        } else {
            checkmate::assertInt(config$max_depth, lower = 1)
        }
        checkmate::assertInt(config$anomalyModelFeatureCount, lower = 1)
        checkmate::assertNumber(config$removeMissingFeatures, lower = 0, upper = 1)
        
        if (!is.null(config$runOrder)) {
            if (!is.data.frame(config$runOrder) && !data.table::is.data.table(config$runOrder)) {
                stop("runOrder must be a data.frame or data.table")
            }
            required_cols = c("Run", "Order")
            missing_cols = setdiff(required_cols, colnames(config$runOrder))
            if (length(missing_cols) > 0) {
                stop("runOrder is missing required columns: ", paste(missing_cols, collapse = ", "))
            }
            if (!is.numeric(config$runOrder$Order)) {
                stop("runOrder$Order must be numeric")
            }
        }
    } else {
        # When anomaly detection is disabled, these parameters are ignored but warn user
        if (length(config$anomalyModelFeatures) > 0) {
            warning("anomalyModelFeatures provided but calculateAnomalyScores=FALSE, ignoring")
        }
        if (!is.null(config$runOrder)) {
            warning("runOrder provided but calculateAnomalyScores=FALSE, ignoring")
        }
    }
    
    # Log file validation
    if (!is.null(config$log_file_path)) {
        checkmate::assertString(config$log_file_path)
        log_dir = dirname(config$log_file_path)
        if (!dir.exists(log_dir)) {
            stop("Log file directory does not exist: ", log_dir)
        }
        if (!file.access(log_dir, mode = 2) == 0) {
            stop("No write permission for log file directory: ", log_dir)
        }
    }
}