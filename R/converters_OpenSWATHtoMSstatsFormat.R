#' Import OpenSWATH files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from OpenSWATH, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param filter_with_mscore TRUE(default) will filter out the features that have greater than mscore_cutoff in m_score column. Those features will be removed.
#' @param mscore_cutoff Cutoff for m_score. Default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#'  
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' @examples 
#' os_raw = system.file("tinytest/raw_data/OpenSWATH/openswath_input.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/annotations/annot_os.csv", 
#'                     package = "MSstats")
#' os_raw = data.table::fread(os_raw) 
#' annot = data.table::fread(annot)
#' 
#' os_imported = OpenSWATHtoMSstatsFormat(os_raw, annot, use_log_file = FALSE)
#' head(os_imported)
#' 
OpenSWATHtoMSstatsFormat = function(
        input, annotation, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
        useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "OpenSWATH", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    m_score_filter = list(score_column = "m_score", 
                          score_threshold = mscore_cutoff, 
                          direction = "smaller", 
                          behavior = "remove", 
                          handle_na = "remove", 
                          fill_value = NA,
                          filter = filter_with_mscore, 
                          drop_column = TRUE)
    
    decoy_filter = list(col_name = "decoy", 
                        filter_symbols = 1, 
                        behavior = "remove",
                        fill_value = NULL,
                        filter = TRUE, 
                        drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        score_filtering = list(ms_filter = m_score_filter),
        exact_filtering = list(decoy = decoy_filter),
        columns_to_fill = c("ProductCharge" = NA, 
                            "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns, 
                                                  fix_missing = "na_to_zero",
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}