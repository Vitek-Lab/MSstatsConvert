#' Import data from SpectroMine
#'
#' @param input data name of SpectroMine PSM output. Read PSM sheet.
#' @param annotation data frame which contains column Run, Fraction, TechRepMixture, Mixture, Channel, BioReplicate, Condition. Refer to the example 'annotation.mine' for the meaning of each column.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with NA and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#'  
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export 
#' @examples
#' raw.mine = data.table::fread(system.file("tinytest/raw_data/SpectroMine/spectromine_input.csv",
#'                                       package = "MSstatsConvert"))
#' annotation.mine = data.table::fread(system.file("tinytest/raw_data/SpectroMine/spectromine_annotation.csv",
#'                                       package = "MSstatsConvert"))
#' head(raw.mine)
#' head(annotation.mine)
#' input.mine <- SpectroMinetoMSstatsTMTFormat(raw.mine, annotation.mine)
#' head(input.mine)
#' 
SpectroMinetoMSstatsTMTFormat <- function(
        input, annotation, filter_with_Qvalue = TRUE, qvalue_cutoff = 0.01,
        useUniquePeptide = TRUE, rmPSM_withfewMea_withinRun = TRUE, 
        rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstatsTMT", "SpectroMine", ...)
    input = MSstatsConvert::MSstatsClean(input)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    pq_filter = list(score_column = "PGQValue", 
                     score_threshold = 0.01, 
                     direction = "smaller",
                     behavior = "fill", 
                     handle_na = "keep", 
                     fill_value = NA,
                     filter = TRUE, 
                     drop_column = TRUE)
    
    qval_filter = list(score_column = "Qvalue", 
                       score_threshold = qvalue_cutoff, 
                       direction = "smaller",
                       behavior = "fill", 
                       handle_na = "keep", 
                       fill_value = NA,
                       filter = filter_with_Qvalue, 
                       drop_column = TRUE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge") 
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        score_filtering = list(pgq = pq_filter, psm_q = qval_filter),
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


