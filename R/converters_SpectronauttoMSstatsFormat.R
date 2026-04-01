#' Import Spectronaut files
#'
#' @param input name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Spectronaut, use annotation=NULL (default). It will use the annotation information from input.
#' @param intensity Intensity column to use. Accepts legacy enum values
#'   \code{'PeakArea'} (default, uses F.PeakArea), \code{'NormalizedPeakArea'}
#'   (uses F.NormalizedPeakArea).  Can also be any raw
#'   Spectronaut column name passed as a string (e.g.
#'   \code{"FG.MS1Quantity"}); the column name is standardized internally.
#'   For protein turnover workflows the recommended default is
#'   \code{"FG.MS1Quantity"}.
#' @param peptideSequenceColumn Name of the Spectronaut column that contains the
#' peptide sequence.  Defaults to \code{"EG.ModifiedSequence"}. The value is 
#' standardized internally (dots and spaces removed) before column lookup.
#' @param heavyLabels Character list identifying the heavy isotope labels as it
#'   appears inside square brackets in the peptide sequence column, e.g.
#'   \code{c("Lys6")} matches peptides containing \code{[Lys6]}.  
#'   \code{c("Lys6", "Arg10")} matches peptides containing either \code{[Lys6]} or \code{[Arg10]}.
#'   Supports any novel label name reported by Spectronaut (e.g. \code{"Leu6"},
#'   \code{"Phe10"}).  When provided, peptides are
#'   classified as heavy (\code{IsotopeLabelType = "H"}), light
#'   (\code{IsotopeLabelType = "L"}), or unlabeled
#'   (\code{IsotopeLabelType = NA}) based on its labeled sequence.  When
#'   \code{NULL} (default) all peptides receive \code{IsotopeLabelType = "L"}.
#'   Useful for protein turnover experiments.
#' @param excludedFromQuantificationFilter Remove rows with F.ExcludedFromQuantification=TRUE Default is TRUE.
#' @param filter_with_Qvalue FALSE(default) will not perform any filtering. TRUE will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param calculateAnomalyScores Default is FALSE. If TRUE, will run anomaly detection model and calculate anomaly scores for each feature. Used downstream to weigh measurements in differential analysis.
#' @param anomalyModelFeatures character vector of quality metric column names to be used as features in the anomaly detection model. List must not be empty if calculateAnomalyScores=TRUE.
#' @param anomalyModelFeatureTemporal character vector of temporal direction corresponding to columns passed to anomalyModelFeatures. Values must be one of: `mean_decrease`, `mean_increase`, `dispersion_increase`, or NULL (to perform no temporal feature engineering). Default is empty vector. If calculateAnomalyScores=TRUE, vector must have as many values as anomalyModelFeatures (even if all NULL).
#' @param removeMissingFeatures Remove features with missing values in more than this fraction of runs. Default is 0.5. Only used if calculateAnomalyScores=TRUE.
#' @param anomalyModelFeatureCount Feature selection for anomaly model. Anomaly detection works on the precursor-level and can be much slower if all features used. We will by default filter to the top-100 highest intensity features. This can be adjusted as necessary. To turn feature-selection off, set this value to a high number (e.g. 10000). Only used if calculateAnomalyScores=TRUE.
#' @param runOrder Temporal order of MS runs. Should be a two column data.table with columns `Run` and `Order`, where `Run` matches the run name output by Spectronaut and `Order` is an integer. Used to engineer the temporal features defined in anomalyModelFeatureTemporal.
#' @param n_trees Number of trees to use in isolation forest when calculateAnomalyScores=TRUE. Default is 100.
#' @param max_depth Max tree depth to use in isolation forest when calculateAnomalyScores=TRUE. Default is "auto" which calculates depth as log2(N) where N is the number of runs. Otherwise must be an integer.
#' @param numberOfCores Number of cores for parallel processing anomaly detection model. When > 1, a logfile named 'MSstats_anomaly_model_progress.log' is created to track progress. Only works for Linux & Mac OS. Default is 1.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
#'                               package = "MSstatsConvert")
#' spectronaut_raw = data.table::fread(spectronaut_raw)
#' spectronaut_imported = SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
#' head(spectronaut_imported)
#' 
SpectronauttoMSstatsFormat = function(
        input, annotation = NULL,
        intensity = 'PeakArea',
        peptideSequenceColumn = "EG.ModifiedSequence",
        heavyLabels = NULL,
        excludedFromQuantificationFilter = TRUE,
        filter_with_Qvalue = FALSE, qvalue_cutoff = 0.01, 
        useUniquePeptide = TRUE, removeFewMeasurements=TRUE,
        removeProtein_with1Feature = FALSE,
        summaryforMultipleRows = max,
        calculateAnomalyScores=FALSE, anomalyModelFeatures=c(),
        anomalyModelFeatureTemporal=c(), removeMissingFeatures=.5,
        anomalyModelFeatureCount=100,
        runOrder=NULL, n_trees=100, max_depth="auto", numberOfCores=1, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, 
        log_file_path = NULL, ...
) {

    validation_config = list(
        input = input,
        annotation = annotation,
        intensity = intensity,
        excludedFromQuantificationFilter = excludedFromQuantificationFilter,
        filter_with_Qvalue = filter_with_Qvalue, 
        qvalue_cutoff = qvalue_cutoff, 
        useUniquePeptide = useUniquePeptide, 
        removeFewMeasurements = removeFewMeasurements,
        removeProtein_with1Feature = removeProtein_with1Feature, 
        summaryforMultipleRows = summaryforMultipleRows, 
        calculateAnomalyScores = calculateAnomalyScores,
        anomalyModelFeatures = anomalyModelFeatures, 
        anomalyModelFeatureTemporal = anomalyModelFeatureTemporal, 
        removeMissingFeatures = removeMissingFeatures,
        anomalyModelFeatureCount = anomalyModelFeatureCount, 
        runOrder = runOrder, 
        n_trees = n_trees, 
        max_depth = max_depth, 
        numberOfCores = numberOfCores,
        use_log_file = use_log_file, 
        append = append, 
        verbose = verbose, 
        log_file_path = log_file_path
    )
    
    .validateMSstatsConverterParameters(validation_config)
    
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    anomalyModelFeatures = .standardizeColnames(anomalyModelFeatures)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Spectronaut", ...)

    input = MSstatsConvert::MSstatsClean(input, intensity = intensity,
                                         calculateAnomalyScores,
                                         anomalyModelFeatures,
                                         peptideSequenceColumn = peptideSequenceColumn,
                                         heavyLabels = heavyLabels)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    pq_filter = list(score_column = "PGQvalue", 
                     score_threshold = 0.01, 
                     direction = "smaller", 
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA_real_,
                     filter = filter_with_Qvalue, 
                     drop_column = TRUE)
    qval_filter = list(score_column = "EGQvalue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller", 
                       behavior = "fill", 
                       handle_na = "keep", 
                       fill_value = NA_real_, 
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    excluded_quant_filter = list(
        col_name    = "FExcludedFromQuantification",
        filter_symbols = TRUE,
        behavior    = "fill",
        fill_value  = NA_real_,
        filter      = excludedFromQuantificationFilter,
        drop_column = TRUE
    )
    
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
    
    fill_isotope_label_type = if (is.null(heavyLabels)) 
        list("IsotopeLabelType" = "L") else list()
    
    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(pgq = pq_filter,
                               psm_q = qval_filter),
        exact_filtering = list(excluded_quant = excluded_quant_filter),
        columns_to_fill = fill_isotope_label_type,
        anomaly_metrics = anomalyModelFeatures)
    input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
    
    input = MSstatsConvert::MSstatsBalancedDesign(
        input, feature_columns, 
        remove_few = removeFewMeasurements,
        anomaly_metrics = anomalyModelFeatures)
    
    if (calculateAnomalyScores){
        input = MSstatsConvert::MSstatsAnomalyScores(
            input, anomalyModelFeatures, anomalyModelFeatureTemporal,
            removeMissingFeatures, anomalyModelFeatureCount, runOrder, n_trees, 
            max_depth, numberOfCores)
    }
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}