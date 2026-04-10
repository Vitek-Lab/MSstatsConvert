#' Generate MSstatsTMT required input format from MaxQuant output
#' 
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param proteinGroups name of 'proteinGroups.txt' data.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mq' for the meaning of each column.
#' @param which.proteinid Use 'Proteins' (default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param rmProt_Only.identified.by.site TRUE will remove proteins with '+' in 'Only.identified.by.site' column from proteinGroups.txt, which was identified only by a modification site. FALSE is the default.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Default is FALSE.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return data.frame of class "MSstatsTMT"
#' 
#' @export
#' 
#' @examples
#' evidence = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_ev.csv",
#'                                       package = "MSstatsConvert"))
#' proteinGroups = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_pg.csv",
#'                                       package = "MSstatsConvert"))
#' annotation.mq = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_annotation.csv",
#'                                       package = "MSstatsConvert"))
#' input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq)
#' head(input.mq)
#' 
MaxQtoMSstatsTMTFormat = function(
        evidence, proteinGroups, annotation, which.proteinid = 'Proteins',
        rmProt_Only.identified.by.site = FALSE, useUniquePeptide = TRUE,
        rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
        summaryforMultipleRows = sum, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    
    input = MSstatsConvert::MSstatsImport(list(evidence = evidence,
                                               protein_groups = proteinGroups), 
                                          "MSstatsTMT", "MaxQuant", ...)
    input = MSstatsConvert::MSstatsClean(
        input,
        protein_id_col = which.proteinid, 
        remove_by_site = rmProt_Only.identified.by.site,
        channel_columns = "Reporterintensitycorrected")
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                                summarize_multiple_psms = summaryforMultipleRows))
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