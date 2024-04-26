#' Import OpenMS files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from OpenMS, which includes feature(peptide ion)-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' 
#' @examples 
#' openms_raw = data.table::fread(system.file("tinytest/raw_data/OpenMS/openms_input.csv", 
#'                                            package = "MSstatsConvert"))
#' openms_imported = OpenMStoMSstatsFormat(openms_raw, use_log_file = FALSE)
#' head(openms_imported)
#' 
OpenMStoMSstatsFormat = function(
        input, annotation = NULL, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "OpenMS", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", 
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}