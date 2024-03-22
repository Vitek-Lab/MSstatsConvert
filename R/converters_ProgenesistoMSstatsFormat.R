#' Import Progenesis files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.
#' @param ... additional parameters to `data.table::fread`.
#'
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek, Ulrich Omasits
#' 
#' @export
#' 
#' @examples 
#' progenesis_raw = system.file("tinytest/raw_data/Progenesis/progenesis_input.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/raw_data/Progenesis/progenesis_annot.csv", 
#'                     package = "MSstatsConvert")
#' progenesis_raw = data.table::fread(progenesis_raw) 
#' annot = data.table::fread(annot)
#' 
#' progenesis_imported = ProgenesistoMSstatsFormat(progenesis_raw, annot,
#'                                                 use_log_file = FALSE)
#' head(progenesis_imported)
#' 
ProgenesistoMSstatsFormat = function(
        input, annotation, useUniquePeptide = TRUE, summaryforMultipleRows = max,
        removeFewMeasurements = TRUE, removeOxidationMpeptides = FALSE, 
        removeProtein_with1Peptide = FALSE, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "Progenesis", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         unique(as.character(annotation$Run)), 
                                         fix_colnames = TRUE)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    oxidation_filter = list(col_name = "PeptideSequence", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        pattern_filtering = list(oxidation = oxidation_filter),
        columns_to_fill = list("FragmentIon" = NA, 
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)
    data.table::setnames(input, "PeptideSequence", "PeptideModifiedSequence",
                         skip_absent = TRUE)
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}