#' Generate MSstatsTMT required input format for OpenMS output
#' @param input MSstatsTMT report from OpenMS
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Default is FALSE.
#' @param summaryforMultiplePSMs sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#'  
#' @return `data.frame` of class `MSstatsTMT`.
#' 
#' @export
#' 
#' @examples
#' raw.om = data.table::fread(system.file("tinytest/raw_data/OpenMSTMT/openmstmt_input.csv",
#'                                       package = "MSstatsConvert"))
#' input.om <- OpenMStoMSstatsTMTFormat(raw.om)
#' head(input.om)
#' 
OpenMStoMSstatsTMTFormat = function(
        input, useUniquePeptide = TRUE, rmPSM_withfewMea_withinRun = TRUE, 
        rmProtein_with1Feature = FALSE, summaryforMultiplePSMs = sum, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstatsTMT", "OpenMS", ...)
    input = MSstatsConvert::MSstatsClean(input)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        NULL, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                                summarize_multiple_psms = summaryforMultiplePSMs)
    )
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  fix_missing = "zero_to_na")
    
    data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the proteinSummarization function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}