#' Import Metamorpheus files
#' 
#' @param input name of Metamorpheus output file, which is tabular format. Use the AllQuantifiedPeaks.tsv file from the Metamorpheus output.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate. 
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Anthony Wu
#' 
#' @export
#' 
#' @examples 
#' input = system.file("tinytest/raw_data/Metamorpheus/AllQuantifiedPeaks.tsv", 
#'                                 package = "MSstatsConvert")
#' input = data.table::fread(input)
#' annot = system.file("tinytest/raw_data/Metamorpheus/Annotation.tsv", 
#'                                 package = "MSstatsConvert")
#' annot = data.table::fread(annot)
#' metamorpheus_imported = MSstatsConvert:::MetamorpheusToMSstatsFormat(input, annotation = annot)
#' head(metamorpheus_imported)
#' 
MetamorpheusToMSstatsFormat = function(
        input, annotation = NULL, useUniquePeptide = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Metamorpheus", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = removeFewMeasurements,
                                summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("FragmentIon" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    return(input)
}