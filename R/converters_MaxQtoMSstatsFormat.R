#' Import MaxQuant files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Raw.file, Condition, BioReplicate, Run, IsotopeLabelType information.
#' @param proteinGroups name of 'proteinGroups.txt' data. It needs to matching protein group ID. If proteinGroups=NULL, use 'Proteins' column in 'evidence.txt'.
#' @param proteinID 'Proteins'(default) or 'Leading.razor.protein' for Protein ID.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#'  
#' @note Warning: MSstats does not support for metabolic labeling or iTRAQ experiments.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 
#' @examples 
#' mq_ev = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_ev.csv",
#'                                       package = "MSstatsConvert"))
#' mq_pg = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_pg.csv",
#'                                       package = "MSstatsConvert"))
#' annot = data.table::fread(system.file("tinytest/raw_data/MaxQuant/annotation.csv",
#'                                       package = "MSstatsConvert"))
#' maxq_imported = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE)
#' head(maxq_imported)
#' 
MaxQtoMSstatsFormat = function(
        evidence, annotation, proteinGroups, proteinID = "Proteins", 
        useUniquePeptide = TRUE, summaryforMultipleRows = max, 
        removeFewMeasurements = TRUE, removeMpeptides = FALSE,
        removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(evidence = evidence, 
                                               protein_groups = proteinGroups), 
                                          type = "MSstats",
                                          tool = "MaxQuant", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         protein_id_col = proteinID, 
                                         remove_by_site = TRUE)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, 
                                                       annotation, 
                                                       "Run" = "Rawfile")
    
    m_filter = list(col_name = "PeptideSequence", 
                    pattern = "M", 
                    filter = removeMpeptides, 
                    drop_column = FALSE)
    
    oxidation_filter = list(col_name = "Modifications", 
                            pattern = "Oxidation", 
                            filter = removeOxidationMpeptides, 
                            drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = useUniquePeptide, 
        remove_single_feature_proteins = removeProtein_with1Peptide,
        pattern_filtering = list(oxidation = oxidation_filter,
                                 m = m_filter),
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
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
    input
}