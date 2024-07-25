#' Generate MSstatsTMT required input format from Protein Prospector output
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input txt report file from Protein Prospector with 
#' `Keep Replicates` option selected.
#' @param annotation data frame which contains column Run, Fraction, 
#' TechRepMixture, Mixture, Channel, BioReplicate, Condition. 
#' 
#' @return data.frame of class "MSstatsTMT"
#' 
#' @export
#' 
#' @examples
#' input = system.file("tinytest/raw_data/ProteinProspector/Prospector_TotalTMT.txt",
#'     package = "MSstatsConvert")
#' input = data.table::fread(input)
#' annot = system.file("tinytest/raw_data/ProteinProspector/Annotation.csv",
#'                                 package = "MSstatsConvert")
#' annot = data.table::fread(annot)
#' output <- ProteinProspectortoMSstatsTMTFormat(input, annot)
#' head(output)
#' 
ProteinProspectortoMSstatsTMTFormat = function(
        input, annotation, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = sum,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, 
        log_file_path = NULL
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstatsTMT", "ProteinProspector")
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows)
        )
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  fix_missing = "zero_to_na")
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the proteinSummarization function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}