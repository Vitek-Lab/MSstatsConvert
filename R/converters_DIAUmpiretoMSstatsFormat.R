#' Import DIA-Umpire files 
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param raw.frag name of FragSummary_date.xls data, which includes feature-level data.
#' @param raw.pep name of PeptideSummary_date.xls data, which includes selected fragments information.
#' @param raw.pro name of ProteinSummary_date.xls data, which includes selected peptides information.
#' @param annotation name of annotation data which includes Condition, BioReplicate, Run information.
#' @param useSelectedFrag TRUE will use the selected fragment for each peptide. 'Selected_fragments' column is required.
#' @param useSelectedPep TRUE will use the selected peptide for each protein. 'Selected_peptides' column is required.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#'
#' @author Meena Choi, Olga Vitek 
#'
#' @export
#' 
#' @examples 
#' diau_frag = system.file("tinytest/raw_data/DIAUmpire/dia_frag.csv", 
#'                              package = "MSstatsConvert")
#' diau_pept = system.file("tinytest/raw_data/DIAUmpire/dia_pept.csv", 
#'                              package = "MSstatsConvert")
#' diau_prot = system.file("tinytest/raw_data/DIAUmpire/dia_prot.csv", 
#'                              package = "MSstatsConvert")
#' annot = system.file("tinytest/raw_data/DIAUmpire/annot_diau.csv", 
#'                     package = "MSstatsConvert")
#' diau_frag = data.table::fread(diau_frag) 
#' diau_pept = data.table::fread(diau_pept) 
#' diau_prot = data.table::fread(diau_prot) 
#' annot = data.table::fread(annot)
#' diau_frag = diau_frag[, lapply(.SD, function(x) if (is.integer(x)) as.numeric(x) else x)]
#' # In case numeric columns are not interpreted correctly
#' 
#' diau_imported = DIAUmpiretoMSstatsFormat(diau_frag, diau_pept, diau_prot, 
#'                                          annot, use_log_file = FALSE)
#' head(diau_imported)
#' 
DIAUmpiretoMSstatsFormat = function(
        raw.frag, raw.pep, raw.pro, annotation, useSelectedFrag = TRUE,
        useSelectedPep = TRUE, removeFewMeasurements = TRUE,
        removeProtein_with1Feature = FALSE, summaryforMultipleRows = max, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(Fragments = raw.frag, 
                                               Peptides = raw.pep, 
                                               Proteins = raw.pro), 
                                          type = "MSstats", 
                                          tool = "DIAUmpire", ...)
    input = MSstatsConvert::MSstatsClean(input, 
                                         use_frag = useSelectedFrag, 
                                         use_pept = useSelectedPep)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "FragmentIon")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation,
        feature_columns,
        remove_shared_peptides = TRUE, 
        remove_single_feature_proteins = removeProtein_with1Feature,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = list("PrecursorCharge" = NA,
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