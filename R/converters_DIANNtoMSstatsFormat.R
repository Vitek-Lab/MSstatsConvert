#' Import Diann files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input name of MSstats input report from Diann, which includes fragment-level data.  
#' Output fragment data with --export-quant flag in DIA-NN 2.0
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' @param MBR True if analysis was done with match between runs
#' @param global_qvalue_cutoff The qvalue cutoff for the Q.Value column, i.e. 
#' the run-specific precursor q-value. Default is 0.01.
#' @param qvalue_cutoff If MBR is false, the qvalue cutoff for the Global.Q.Value 
#' column, i.e. global precursor q-value.  If MBR is true, the qvalue cutoff for the 
#' Lib.Q.Value column, i.e. the q-value for the library created after the first MBR pass.
#' Default is 0.01.
#' @param pg_qvalue_cutoff If MBR is false, the qvalue cutoff for the Global.PG.Q.Value 
#' column, i.e. the global q-value for the protein group.  If MBR is true, the
#' qvalue cutoff for the Lib.PG.Q.Value column, i.e. the protein group q-value for 
#' the library created after the first MBR pass. Default is 0.01.
#' @param useUniquePeptide should unique pepties be removed
#' @param removeFewMeasurements should proteins with few measurements be removed
#' @param removeOxidationMpeptides should peptides with oxidation be removed
#' @param removeProtein_with1Feature should proteins with a single feature be removed
#' @param quantificationColumn Use 'FragmentQuantCorrected'(default) column for quantified intensities for DIANN 1.8.x.
#' Use 'FragmentQuantRaw' for quantified intensities for DIANN 1.9.x. 
#' Use 'auto' for quantified intensities for DIANN 2.x where each fragment intensity is a separate column, e.g. Fr0Quantity.
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
#' 
#' # For DIANN 2.0, set quantificationColumn = 'auto'
#' input_file_path_2_0 = system.file("tinytest/raw_data/DIANN/diann_2.0.parquet", 
#'                                 package="MSstatsConvert")
#' annotation_file_path_2_0 = system.file("tinytest/raw_data/DIANN/annotation_diann_2.0.csv", 
#'                                 package = "MSstatsConvert")
#' input_2_0 = arrow::read_parquet(input_file_path_2_0)
#' annot_2_0 = data.table::fread(annotation_file_path_2_0)
#' output_2_0 = DIANNtoMSstatsFormat(input_2_0, annotation = annot_2_0, MBR = FALSE, 
#'                                 use_log_file = FALSE, quantificationColumn = 'auto')
#' head(output_2_0)
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
                                MBR = TRUE, 
                                quantificationColumn = "FragmentQuantCorrected",
                                ...) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "DIANN")
    input = MSstatsConvert::MSstatsClean(input, MBR = MBR, 
                                         quantificationColumn = quantificationColumn,
                                         global_qvalue_cutoff = global_qvalue_cutoff,
                                         qvalue_cutoff = qvalue_cutoff, 
                                         pg_qvalue_cutoff = pg_qvalue_cutoff)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)
    
    decoy_filter = list(col_name = "ProteinName",
                        pattern = c("DECOY", "Decoys"),
                        filter = T, 
                        drop_column = FALSE)
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "\\(UniMod\\:35\\)",
                            filter = removeOxidationMpeptides,
                            drop_column = FALSE)
    
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
    input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]
    
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