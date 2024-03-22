#' Import FragPipe files
#' 
#' @param input name of FragPipe msstats.csv export. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity are required.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Devon Kohler
#' 
#' @export
#' 
#' @examples 
#' fragpipe_raw = system.file("tinytest/raw_data/FragPipe/fragpipe_input.csv",
#'                               package = "MSstatsConvert")
#' fragpipe_raw = data.table::fread(fragpipe_raw)
#' fragpipe_imported = FragPipetoMSstatsFormat(fragpipe_raw, use_log_file = FALSE)
#' head(fragpipe_imported)
#' 
FragPipetoMSstatsFormat = function(
        input, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "FragPipe", ...)
    input = MSstatsConvert::getInputFile(input, "input")
    
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, NULL)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}