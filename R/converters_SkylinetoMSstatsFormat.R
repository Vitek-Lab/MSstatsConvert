#' Import Skyline files
#'
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeled 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' skyline_raw = system.file("tinytest/raw_data/Skyline/skyline_input.csv",
#'                           package = "MSstatsConvert")
#' skyline_raw = data.table::fread(skyline_raw)
#' skyline_imported = SkylinetoMSstatsFormat(skyline_raw)
#' head(skyline_imported)
#' 
SkylinetoMSstatsFormat = function(
        input, annotation = NULL, removeiRT = TRUE, filter_with_Qvalue = TRUE,
        qvalue_cutoff = 0.01, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Skyline", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, 
                                                       annotation, 
                                                       Run = "FileName")
    
    decoy_filter = list(col_name = "ProteinName",
                        pattern = c("DECOY", "Decoys"),
                        filter = TRUE, 
                        drop_column = FALSE)
    
    irt_filter = list(col_name = "StandardType", 
                      filter_symbols = "iRT",
                      filter = removeiRT, 
                      behavior = "remove",
                      fill_value = NULL,
                      drop_column = FALSE)
    
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\+16", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    truncated_filter = list(col_name = "Truncated", 
                            filter_symbols = "TRUE",
                            behavior = "fill",
                            fill_value = NA_real_,
                            filter = TRUE, 
                            drop_column = TRUE)
    
    qval_filter = list(score_column = "DetectionQValue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller",
                       behavior = "fill", 
                       fill_value = 0, 
                       handle_na = "keep",
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    
    feature_columns = c("IsotopeLabelType", "PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        score_filtering = list(qval = qval_filter),
        pattern_filtering = list(decoy = decoy_filter, 
                                 oxidation = oxidation_filter),
        exact_filtering = list(irt = irt_filter,
                               truncated = truncated_filter),
        aggregate_isotopic = TRUE,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = sum))
    input = MSstatsBalancedDesign(input, c("PeptideSequence", "PrecursorCharge", 
                                           "FragmentIon", "ProductCharge"),
                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}