#' Convert Philosopher (Fragpipe) output to MSstatsTMT format.
#' 
#' @param input data.frame of `msstats.csv` file produced by Philosopher
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column. Channel column should be 
#' consistent with the channel columns (Ignore the prefix "Channel ") in msstats.csv file. Run column should be consistent with the Spectrum.File columns in msstats.csv file.
#' @param protein_id_col Use 'Protein'(default) column for protein name. 
#' 'Master.Protein.Accessions' can be used instead to get the protein ID with single protein.
#' @param peptide_id_col Use 'Peptide.Sequence'(default) column for peptide sequence.
#'  'Modified.Peptide.Sequence' can be used instead to get the modified peptide sequence.
#' @param Purity_cutoff Cutoff for purity. Default is 0.6
#' @param PeptideProphet_prob_cutoff Cutoff for the peptide identification probability. Default is 0.7. 
#' The probability is confidence score determined by PeptideProphet and higher values indicate greater confidence.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmPeptide_OxidationM TRUE (default) will remove the peptides including oxidation (M) sequence.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Default is FALSE.
#' @param ... additional parameters to `data.table::fread`.
#' @inheritParams .sharedParametersAmongConverters
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @examples
#' input_file_path = system.file("tinytest/raw_data/Philosopher/msstats.csv", 
#'                      package = "MSstatsConvert")
#' annotation_file_path = system.file("tinytest/raw_data/Philosopher/MSstatsTMT_annotation.csv", 
#'                     package = "MSstatsConvert")
#' input = data.table::fread(input_file_path)
#' annotation = data.table::fread(annotation_file_path)
#' msstats_format = PhilosophertoMSstatsTMTFormat(input, annotation)
#' head(msstats_format)
#' 
#' @export
PhilosophertoMSstatsTMTFormat = function(
        input, annotation, protein_id_col = "Protein", 
        peptide_id_col = "Peptide.Sequence", Purity_cutoff = 0.6, 
        PeptideProphet_prob_cutoff = 0.7, useUniquePeptide = TRUE, 
        rmPSM_withfewMea_withinRun = TRUE, rmPeptide_OxidationM = TRUE, 
        rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum, 
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path, 
                                        base = "MSstatsTMT_converter_log_")
    checkmate::assertTRUE(!is.null(input))
    
    input[["Is.Unique"]] = as.logical(input[["Is.Unique"]])
    mixture_files = list(Mixture1=data.table(input))
    
    input = MSstatsConvert::MSstatsImport(c(mixture_files,
                                            list(annotation = annotation)), 
                                          type = "MSstatsTMT",
                                          tool = "Philosopher")
    channels = unique(annotation[["Channel"]])
    input = MSstatsConvert::MSstatsClean(input, protein_id_col, peptide_id_col,
                                         channels, useUniquePeptide)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    purity_filter = list(score_column = "Purity", 
                         score_threshold = Purity_cutoff, 
                         direction = "greater",
                         behavior = "remove", 
                         handle_na = "keep", 
                         fill_value = NULL,
                         filter = TRUE, 
                         drop_column = FALSE)
    probability_filter = list(score_column = "PeptideProphetProbability", 
                              score_threshold = PeptideProphet_prob_cutoff, 
                              direction = "greater",
                              behavior = "remove", 
                              handle_na = "keep", 
                              fill_value = NULL,
                              filter = TRUE, 
                              drop_column = FALSE)
    oxidation_filter = list(col_name = "PeptideSequence",
                            pattern = "Oxidation", 
                            filter = rmPeptide_OxidationM, 
                            drop_column = FALSE)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge") 
    input = MSstatsPreprocess(
        input, 
        annotation, 
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        score_filtering = list(pur = purity_filter, prob = probability_filter),
        pattern_filtering = list(ox = oxidation_filter),
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


#' Convert Philosopher parameters to consistent format
#' @inheritParams PhilosophertoMSstatsTMTFormat
#' @param path character. Path to a file or directory containing msstats.csv output(s) from Philosopher. Used when \code{input} is NULL.
#' @param folder logical. If TRUE, \code{path} is treated as a directory and all msstats files within it are read. If FALSE, \code{path} is treated as a single file path.
#' @keywords internal
.getPhilosopherInput = function(input, path, folder) {
    if (!is.null(input)) {
        mixture_files = input
    } else {
        if (folder) {
            mixture_files = list.files(path, pattern = "msstats", 
                                       full.names = TRUE)
        } else {
            mixture_files = path
        }
    }
    mixture_files = as.list(mixture_files)
    names(mixture_files) = paste0("Mixture", seq_along(mixture_files))    
    mixture_files
}
