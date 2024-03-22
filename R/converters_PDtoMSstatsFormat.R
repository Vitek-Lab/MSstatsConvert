#' Import Proteome Discoverer files
#' 
#' @inheritParams .sharedParametersAmongConverters
#' @param input PD report or a path to it.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, 
#' Run information. 'Run' will be matched with 'Spectrum.File'.
#' @param useNumProteinsColumn TRUE removes peptides which have more than 1 in # Proteins column of PD output.
#' @param which.quantification Use 'Precursor.Area'(default) column for quantified intensities. 'Intensity' or 'Area' can be used instead.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead.
#' @param which.sequence Use 'Sequence'(default) column for peptide sequence. 'Annotated.Sequence' can be used instead.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return data.frame in the MSstats required format.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
#' @examples 
#' 
#' pd_raw = system.file("tinytest/raw_data/PD/pd_input.csv", 
#'                      package = "MSstatsConvert")
#' annot = system.file("tinytest/annotations/annot_pd.csv", package = "MSstats")
#' pd_raw = data.table::fread(pd_raw)
#' annot = data.table::fread(annot)
#' 
#' pd_imported = PDtoMSstatsFormat(pd_raw, annot, use_log_file = FALSE)
#' head(pd_imported)
#' 
PDtoMSstatsFormat = function(
        input, annotation, useNumProteinsColumn = FALSE, useUniquePeptide = TRUE,
        summaryforMultipleRows = max, removeFewMeasurements = TRUE,
        removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
        which.quantification = 'Precursor.Area', 
        which.proteinid = 'Protein.Group.Accessions', which.sequence = 'Sequence',
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose, 
                                        log_file_path)
    
    input = MSstatsConvert::MSstatsImport(list(input = input), 
                                          "MSstats", "ProteomeDiscoverer", ...)
    input = MSstatsConvert::MSstatsClean(
        input, 
        quantification_column = which.quantification, 
        protein_id_column = which.proteinid,
        sequence_column = which.sequence, 
        remove_shared = useNumProteinsColumn)
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