#' Convert Proteome Discoverer output to MSstatsTMT format.
#' 
#' @param input PD report or a path to it.
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead to get the protein name with single protein.
#' @param useNumProteinsColumn logical, TRUE (default) removes shared peptides by information of # Proteins column in PSM sheet.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Default is FALSE.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export
#' 
#' @examples
#' raw.pd = data.table::fread(system.file("tinytest/raw_data/PDTMT/pdtmt_input.csv",
#'                                       package = "MSstatsConvert"))
#' annotation.pd = data.table::fread(system.file("tinytest/raw_data/PDTMT/pd_annotation.csv",
#'                                       package = "MSstatsConvert"))
#' head(raw.pd)
#' head(annotation.pd)
#' input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)
#' head(input.pd)
#' 
PDtoMSstatsTMTFormat <- function(
        input, annotation, which.proteinid = 'Protein.Accessions', 
        useNumProteinsColumn = TRUE, useUniquePeptide = TRUE, 
        rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
        summaryforMultipleRows = sum, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    
    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstatsTMT", "ProteomeDiscoverer", ...)
    input = MSstatsConvert::MSstatsClean(
        input, 
        protein_id_column = which.proteinid,
        remove_shared = useNumProteinsColumn,
        remove_protein_groups = useNumProteinsColumn)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation, 
        feature_columns = feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                                summarize_multiple_psms = summaryforMultipleRows)
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