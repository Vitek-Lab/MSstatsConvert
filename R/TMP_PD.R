#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' @importFrom MSstatsConvert MSstatsImport MSstatsClean MSstatsPreprocess 
#' MSstatsBalancedDesign MSstatsMakeAnnotation MSstatsSaveSessionInfo
#' MSstatsLogsSettings
#' 
#' @param fewMeasurements 'remove'(default) will remove the features that have 1 or 2 measurements across runs.
#' @param useUniquePeptide TRUE (default) removes peptides that are assigned for more than one proteins. 
#' We assume to use unique peptide for each protein.
#' @param summaryforMultipleRows max(default) or sum - when there are multiple measurements for certain feature and certain run, use highest or sum of multiple intensities.
#' @param removeProtein_with1Feature TRUE will remove the proteins which have only 1 feature, which is the combination of peptide, precursor charge, fragment and charge. FALSE is default.
#' @param removeProtein_with1Peptide TRUE will remove the proteins which have only 1 peptide and charge. FALSE is default.
#' @param removeOxidationMpeptides TRUE will remove the peptides including 'oxidation (M)' in modification. FALSE is default.
#' @param removeMpeptides TRUE will remove the peptides including 'M' sequence. FALSE is default.
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be added
#' to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be printed
#' to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' @param sesion_info_path character. Optional path to a file to which session 
#' information will be saved.
#' 
#' @keywords internal
#' 

.documentFunction = function(fewMeasurements, 
                             useUniquePeptide,
                             summaryforMultipleRows, 
                             removeProtein_with1Feature,
                             removeProtein_with1Protein,
                             removeOxidationMpeptides,
                             removeMpeptides) {
    NULL
}


standard_columns_tmt = c("ProteinName", "PeptideSequence", "Charge", "PSM", 
                         "Mixture", "TechRepMixture", "Run", "Channel", 
                         "BioReplicate", "Condition", "Intensity" )


#' Convert Proteome Discoverer output to MSstatsTMT format.
#' 
#' @param input PD report or a path to it.
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param useNumProteinsColumn logical, if TRUE, shared peptides will be removed.
#' @param useUniquePeptide lgl, if TRUE (default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withfewMea_withinRun TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum (default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @param ... additional parameters to `data.table::fread`.
#' 
#' @return `data.frame` of class `MSstatsTMT`
#' 
#' @export
#' 
PDtoMSstatsTMTFormat <- function(
    input, annotation, which.proteinid = 'Protein.Accessions', 
    useNumProteinsColumn = TRUE, useUniquePeptide = TRUE, 
    rmPSM_withfewMea_withinRun = TRUE, rmProtein_with1Feature = FALSE, 
    summaryforMultipleRows = sum, 
    use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
    session_info_path = NULL, ...
) {
    MSstatsLogsSettings(use_log_file, append, verbose, log_file_path,
                        base = "MSstatsTMT_log_")
    MSstatsSaveSessionInfo(session_info_path, append = TRUE)
    
    input = MSstatsImport(list(input = input),
                          "MSstatsTMT", "ProteomeDiscoverer", ...)
    input = MSstatsClean(input, 
                         protein_id_column = which.proteinid,
                         remove_shared = useNumProteinsColumn,
                         remove_protein_groups = useNumProteinsColumn)
    annotation = MSstatsMakeAnnotation(input, annotation)
    
    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsPreprocess(
        input,
        annotation, 
        feature_columns = c("PeptideSequence", "PrecursorCharge"),
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = rmProtein_with1Feature,
        feature_cleaning = list(remove_features_with_few_measurements = rmPSM_withfewMea_withinRun,
                                summarize_multiple_psms = summaryforMultipleRows)
    )
    input = MSstatsBalancedDesign(input, feature_columns)
    data.table::setnames(input, "PrecursorCharge", "Charge", skip_absent = TRUE)
    input[, intersect(standard_columns_tmt, colnames(input))]
}
