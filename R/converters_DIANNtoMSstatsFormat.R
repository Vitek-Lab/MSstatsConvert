#' Import Diann files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from Diann, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' @param MBR True if analysis was done with match between runs
#' @param global_qvalue_cutoff The global qvalue cutoff
#' @param qvalue_cutoff local qvalue cutoff for library
#' @param pg_qvalue_cutoff local qvalue cutoff for protein groups Run should be the same as filename.
#' @param useUniquePeptide should unique pepties be removed
#' @param removeFewMeasurements should proteins with few measurements be removed
#' @param removeOxidationMpeptides should peptides with oxidation be removed
#' @param removeProtein_with1Feature should proteins with a single feature be removed
#' @param ... additional parameters to `data.table::fread`.
#'  
#' @return data.frame in the MSstats required format.
#' 
#' @author Elijah Willie
#' 
#' @export
#' 
#' @examples 
#' input_file_path = system.file("tinytest/raw_data/DIANN/diann_input.tsv", 
#'                                 package="MSstatsConvert")
#' annotation_file_path = system.file("tinytest/raw_data/DIANN/annotation.csv", 
#'                                 package = "MSstatsConvert")
#' input = data.table::fread(input_file_path)
#' annot = data.table::fread(annotation_file_path)
#' output = DIANNtoMSstatsFormat(input, annotation = annot, MBR = FALSE, 
#'                                 use_log_file = FALSE)
#' head(output)
DIANNtoMSstatsFormat = function(input, annotation = NULL,
                                global_qvalue_cutoff = 0.01,
                                qvalue_cutoff = 0.01, 
                                pg_qvalue_cutoff = 0.01,
                                useUniquePeptide = TRUE, 
                                removeFewMeasurements = TRUE,
                                removeOxidationMpeptides = TRUE, 
                                removeProtein_with1Feature = TRUE,
                                use_log_file = TRUE, append = FALSE, 
                                verbose = TRUE, log_file_path = NULL,
                                MBR = TRUE,...) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "DIANN")
    input = MSstatsConvert::MSstatsClean(input, MBR = MBR)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    decoy_filter = list(col_name = "ProteinName",
                        pattern = c("DECOY", "Decoys"),
                        filter = T, 
                        drop_column = FALSE)
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\(UniMod\\:35\\)",
                            filter = removeOxidationMpeptides,
                            drop_column = FALSE)
    
    msg = paste0('** Filtering on Global Q Value < ', global_qvalue_cutoff)
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    
    input = input[DetectionQValue < global_qvalue_cutoff, ]
    if (MBR) {
        msg = '** MBR was used to analyze the data. Now setting names and filtering'
        msg_1_mbr = paste0('-- LibPGQValue < ', pg_qvalue_cutoff)
        msg_2_mbr = paste0('-- LibQValue < ', qvalue_cutoff)
        input = input[LibPGQValue < pg_qvalue_cutoff, ]
        input = input[LibQValue < qvalue_cutoff, ]
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg_1_mbr)
        getOption("MSstatsMsg")("INFO", msg_1_mbr)
        getOption("MSstatsLog")("INFO", msg_2_mbr)
        getOption("MSstatsMsg")("INFO", msg_2_mbr)
        # getOption("MSstatsLog")("INFO", "\n")
    } else{
        msg = '** MBR was not used to analyze the data. Now setting names and filtering'
        msg_1 = paste0('-- Filtering on GlobalPGQValue < ', pg_qvalue_cutoff)
        msg_2 = paste0('-- Filtering on GlobalQValue < ', qvalue_cutoff)
        input = input[GlobalPGQValue < pg_qvalue_cutoff, ]
        input = input[GlobalQValue < qvalue_cutoff, ]
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg_1)
        getOption("MSstatsMsg")("INFO", msg_1)
        getOption("MSstatsLog")("INFO", msg_2)
        getOption("MSstatsMsg")("INFO", msg_2)
        # getOption("MSstatsLog")("INFO", "\n")
    }
    
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Feature,
        exact_filtering = NULL,
        pattern_filtering = list(decoy = decoy_filter, 
                                 oxidation = oxidation_filter),
        aggregate_isotopic = FALSE,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = max),
        columns_to_fill = list(Fraction = 1,
                               IsotopeLabelType = "Light"))
    
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns, 
                                                  fill_incomplete = FALSE,
                                                  handle_fractions = FALSE,
                                                  remove_few = removeFewMeasurements
    )
    
    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}