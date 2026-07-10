#' Import Diann files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from Diann, which includes fragment-level data.  
#' Output fragment data with --export-quant flag in DIA-NN 2.0
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' @param MBR True if analysis was done with match between runs
#' @param global_qvalue_cutoff The qvalue cutoff for the Q.Value column, i.e. 
#' the run-specific precursor q-value. Default is 0.01.
#' @param qvalue_cutoff If MBR is false, the qvalue cutoff for the Global.Q.Value 
#' column, i.e. global precursor q-value.  If MBR is true, the qvalue cutoff for the 
#' Lib.Q.Value column, i.e. the q-value for the library created after the first MBR pass.
#' Default is 0.01.
#' @param pg_qvalue_cutoff If MBR is false, the qvalue cutoff for the Global.PG.Q.Value 
#' column, i.e. the global q-value for the protein group.  If MBR is true, the
#' qvalue cutoff for the Lib.PG.Q.Value column, i.e. the protein group q-value for 
#' the library created after the first MBR pass. Default is 0.01.
#' @param useUniquePeptide should unique peptides be removed
#' @param removeFewMeasurements should proteins with few measurements be removed
#' @param removeOxidationMpeptides should peptides with oxidation be removed
#' @param removeProtein_with1Feature should proteins with a single feature be removed
#' @param labeledAminoAcids Character vector of single-letter amino acid codes
#' that carry the SILAC label in protein turnover experiments, e.g.
#' \code{c("K")} or \code{c("K", "R")}.  Supplying this vector opts in to
#' protein-turnover mode; the exact amino acids determine behaviour only in the
#' \code{ModifiedSequence}-parsing path described below.
#'
#' \strong{Channel-based path} (DIA-NN 2.x exports that include a
#' \code{Channel} column): when \code{labeledAminoAcids} is non-\code{NULL}
#' \emph{and} the input contains a \code{Channel} column, \code{Channel} values
#' are mapped directly to \code{IsotopeLabelType} (\code{"H"} → \code{"H"},
#' \code{"L"} → \code{"L"}, anything else → \code{NA}).  The amino acid codes
#' in \code{labeledAminoAcids} are \strong{not} used to validate or filter
#' \code{ModifiedSequence} in this path.
#'
#' \strong{ModifiedSequence-parsing path} (DIA-NN 1.x exports without a
#' \code{Channel} column): when \code{labeledAminoAcids} is non-\code{NULL}
#' and no \code{Channel} column is present, each \code{ModifiedSequence} is 
#' scanned for isotope-labeled amino acids, which appear in parentheses 
#' immediately after the labeled residue, 
#' in the form \code{(<label>-<aminoAcid>-H)} for
#' the heavy label or \code{(<label>-<aminoAcid>-L)} for the light label.
#' For example, \code{K(label-K-H)} marks a heavy-labeled lysine (\code{K}).
#' Here \code{<aminoAcid>} is one of the single-letter codes in
#' \code{labeledAminoAcids}, and \code{<label>} is the label name (e.g.
#' \code{SILAC} or \code{label}). Sequences with no such parenthetical
#' tags are assigned \code{IsotopeLabelType = NA}. Once classified, the
#' parenthetical annotation is stripped out of \code{PeptideSequence},
#' leaving the plain amino acid sequence.
#'
#' When \code{NULL} (default), protein-turnover mode is disabled and all
#' peptides receive \code{IsotopeLabelType = "Light"}.
#' @param quantificationColumn Use 'FragmentQuantCorrected'(default) column for quantified intensities for DIANN 1.8.x.
#' Use 'FragmentQuantRaw' for quantified intensities for DIANN 1.9.x.
#' Use 'auto' for quantified intensities for DIANN 2.x where each fragment intensity is a separate column, e.g. Fr0Quantity.
#' @param calculateAnomalyScores Default is FALSE. If TRUE, will run anomaly detection model and calculate anomaly scores for each feature. Used downstream to weigh measurements in differential analysis.
#' @param anomalyModelFeatures character vector of quality metric column names to be used as features in the anomaly detection model. List must not be empty if calculateAnomalyScores=TRUE.
#' @param anomalyModelFeatureTemporal character vector of temporal direction corresponding to columns passed to anomalyModelFeatures. Values must be one of: `mean_decrease`, `mean_increase`, `dispersion_increase`, or NULL (to perform no temporal feature engineering). Default is empty vector. If calculateAnomalyScores=TRUE, vector must have as many values as anomalyModelFeatures (even if all NULL).
#' @param removeMissingFeatures Remove features with missing values in more than this fraction of runs. Default is 0.5. Only used if calculateAnomalyScores=TRUE.
#' @param anomalyModelFeatureCount Feature selection for anomaly model. Anomaly detection works on the precursor-level and can be much slower if all features used. We will by default filter to the top-100 highest intensity features. This can be adjusted as necessary. To turn feature-selection off, set this value to a high number (e.g. 10000). Only used if calculateAnomalyScores=TRUE.
#' @param runOrder Temporal order of MS runs. Should be a two column data.table with columns `Run` and `Order`, where `Run` matches the run name output by DIA-NN and `Order` is an integer. Used to engineer the temporal features defined in anomalyModelFeatureTemporal.
#' @param n_trees Number of trees to use in isolation forest when calculateAnomalyScores=TRUE. Default is 100.
#' @param max_depth Max tree depth to use in isolation forest when calculateAnomalyScores=TRUE. Default is "auto" which calculates depth as log2(N) where N is the number of runs. Otherwise must be an integer.
#' @param numberOfCores Number of cores for parallel processing anomaly detection model. When > 1, a logfile named 'MSstats_anomaly_model_progress.log' is created to track progress. Only works for Linux & Mac OS. Default is 1.
#' @param ... additional parameters to `data.table::fread`.
#'  
#' @return data.frame in the MSstats required format.
#' 
#' @author Elijah Willie
#' 
#' @export
#' 
#' @examples 
#' input_file_path = system.file("tinytest/raw_data/DIANN/diann_input.tsv", 
#'                                 package="MSstatsConvert")
#' annotation_file_path = system.file("tinytest/raw_data/DIANN/annotation.csv", 
#'                                 package = "MSstatsConvert")
#' input = data.table::fread(input_file_path)
#' annot = data.table::fread(annotation_file_path)
#' output = DIANNtoMSstatsFormat(input, annotation = annot, MBR = FALSE, 
#'                                 use_log_file = FALSE)
#' head(output)
#' 
#' # For DIANN 2.0, set quantificationColumn = 'auto'
#' input_file_path_2_0 = system.file("tinytest/raw_data/DIANN/diann_2.0.parquet", 
#'                                 package="MSstatsConvert")
#' annotation_file_path_2_0 = system.file("tinytest/raw_data/DIANN/annotation_diann_2.0.csv", 
#'                                 package = "MSstatsConvert")
#' input_2_0 = arrow::read_parquet(input_file_path_2_0)
#' annot_2_0 = data.table::fread(annotation_file_path_2_0)
#' output_2_0 = DIANNtoMSstatsFormat(input_2_0, annotation = annot_2_0, MBR = FALSE, 
#'                                 use_log_file = FALSE, quantificationColumn = 'auto')
#' head(output_2_0)
DIANNtoMSstatsFormat = function(
    input, annotation = NULL,
    global_qvalue_cutoff = 0.01,
    qvalue_cutoff = 0.01,
    pg_qvalue_cutoff = 0.01,
    useUniquePeptide = TRUE,
    removeFewMeasurements = TRUE,
    removeOxidationMpeptides = TRUE,
    removeProtein_with1Feature = TRUE,
    MBR = TRUE,
    labeledAminoAcids = NULL,
    quantificationColumn = "FragmentQuantCorrected",
    calculateAnomalyScores=FALSE, anomalyModelFeatures=c(),
    anomalyModelFeatureTemporal=c(), removeMissingFeatures=.5,
    anomalyModelFeatureCount=100,
    runOrder=NULL, n_trees=100, max_depth="auto", numberOfCores=1,
    use_log_file = TRUE, append = FALSE,
    verbose = TRUE, log_file_path = NULL,
    ...) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    anomalyModelFeatures = .standardizeColnames(anomalyModelFeatures)

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "DIANN")
   
    input = MSstatsConvert::MSstatsClean(input, MBR = MBR,
                                         quantificationColumn = quantificationColumn,
                                         global_qvalue_cutoff = global_qvalue_cutoff,
                                         qvalue_cutoff = qvalue_cutoff,
                                         pg_qvalue_cutoff = pg_qvalue_cutoff,
                                         calculateAnomalyScores = calculateAnomalyScores,
                                         anomalyModelFeatures = anomalyModelFeatures,
                                         labeledAminoAcids = labeledAminoAcids)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

    decoy_filter = list(col_name = "ProteinName",
                        pattern = c("DECOY", "Decoys"),
                        filter = T,
                        drop_column = FALSE)
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\(UniMod\\:35\\)",
                            filter = removeOxidationMpeptides,
                            drop_column = FALSE)

    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
    preprocess_feature_columns = if ("IsotopeLabelType" %in% colnames(input))
        c(feature_columns, "IsotopeLabelType") else feature_columns
    fill_isotope_label_type = if ("IsotopeLabelType" %in% colnames(input))
        list() else list("IsotopeLabelType" = "Light")

    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        preprocess_feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        exact_filtering = NULL,
        pattern_filtering = list(decoy = decoy_filter,
                                 oxidation = oxidation_filter),
        aggregate_isotopic = FALSE,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = max),
        columns_to_fill = c(list(Fraction = 1), fill_isotope_label_type),
        anomaly_metrics = anomalyModelFeatures)
    input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
    # browser()
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns, 
                                                  fill_incomplete = TRUE,
                                                  handle_fractions = FALSE,
                                                  remove_few = removeFewMeasurements,
                                                  anomaly_metrics = anomalyModelFeatures
    )
    # browser()
    if (calculateAnomalyScores){
        input = MSstatsConvert::MSstatsAnomalyScores(
            input, anomalyModelFeatures, anomalyModelFeatureTemporal,
            removeMissingFeatures, anomalyModelFeatureCount, runOrder, n_trees, 
            max_depth, numberOfCores)
    }
    # browser()
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}
