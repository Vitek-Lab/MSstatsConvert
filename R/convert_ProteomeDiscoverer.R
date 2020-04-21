#' Import Proteome Discoverer files
#' 
#' @inheritParams .documentFunction
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, 
#' Run information. 'Run' will be matched with 'Spectrum.File'.
#' @param useNumProteinsColumn TRUE removes peptides which have more than 1 in # Proteins column of PD output.
#' @param which.quantification Use 'Precursor.Area'(default) column for quantified intensities. 'Intensity' or 'Area' can be used instead.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead.
#' @param which.sequence Use 'Sequence'(default) column for peptide sequence. 'Annotated.Sequence' can be used instead.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 

.PDFromDFs = function(
    input, annotation, useNumProteinsColumn = FALSE, useUniquePeptide = TRUE,
    summaryforMultipleRows = max, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    which.quantification = 'Precursor.Area', 
    which.proteinid = 'Protein.Group.Accessions', which.sequence = 'Sequence',
    use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
    
    fewMeasurements = .isLegalValue(fewMeasurements, 
                                    legal_values = c("remove", "keep"))
    
    annotation = .makeAnnotation(annotation, 
                                 c("Run" = "Run", "Condition" = "Condition", 
                                   "BioReplicate" = "BioReplicate"))
    input = .cleanRawPD(input, which.proteinid, which.quantification, 
                        which.sequence, useNumProteinsColumn)
    feature_cols = c("PeptideModifiedSequence", "PrecursorCharge", "FragmentIon")
    input = .handleOxidationPeptides(input, "PeptideModifiedSequence", 
                                     "Oxidation", removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide,
                                  peptide_column = "PeptideModifiedSequence")
    input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                            fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide)
    input = merge(input, annotation, by = "Run")
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 "IsotopeLabelType" = "L"))
    input
}


.PDFromFiles = function(
    input, annotation, useNumProteinsColumn = FALSE, useUniquePeptide = TRUE,
    summaryforMultipleRows = max, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    which.quantification = 'Precursor.Area', 
    which.proteinid = 'Protein.Group.Accessions', which.sequence = 'Sequence',
    use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
    
}

setGeneric("PDtoMSstatsFormat",
           function(input, annotation, ...) {
               standardGeneric("PDtoMSstatsFormat")
           })

setMethod("PDtoMSstatsFormat",
          signature = rep("character", 2),
          function(input, annotation, ...)
              .PDFromFiles(input, annotation, ...))


setMethod("PDtoMSstatsFormat",
          signature = rep("data.frame", 2),
          function(input, annotation, ...) {
              .PDFromDFs(input, annotation, ...)
          })



.cleanRawPD = function(pd_input, quantification_column, proteinID_column,
                       sequence_column, filter_num_col) {
    colnames(pd_input) = .standardizeColnames(pd_input)
    which.quantification = .findAvailable(c("Intensity", "Area"),
                                          colnames(pd_input),
                                          "Precursor.Area")
    which.quantification = .isLegalValue(which.quantification, 
                                         legal_values = c("Intensity", "Area", "Precursor.Area",
                                                          "Precursor.Abundance"),
                                         message = "Please select a column to be used for quantified intensities among four options: ")
    
    which.proteinid = .findAvailable(c("Protein.Accessions", 
                                       "Master.Protein.Accessions"),
                                     colnames(pd_input),
                                     "Protein.Group.Accessions")
    which.proteinid = .isLegalValue(which.proteinid, 
                                    legal_values = c("ProteinAccessions", 
                                                     "Master.Protein.Accessions",
                                                     "Protein.Group.Accessions"),
                                    message = "Please select a column to be used as protein IDs among three options: ")
    which.sequence = .findAvailable("Annotated.Sequence", colnames(pd_input), 
                                    "Sequence")
    which.sequence = .isLegalValue(which.sequence, 
                                   legal_values = c("Annotated.Sequence", "Sequence"),
                                   message = "Please select peptide sequence column between two options: ")
    
    if(filter_num_col) {
        pd_input = pd_input[pd_input[["X..Proteins"]] == '1', ]
    }
    pd_cols = c(proteinID_column, "X..Proteins", sequence_column, 
                "Modifications", "Charge", "Spectrum.File", quantification_column)
    if (any(is.element(colnames(pd_input), "Fraction"))) {
        pd_cols = c(pd_cols, "Fraction")
    }
    input = data.table::as.data.table(pd_input[, pd_cols])
    colnames(input) = .updateColnames(
        input,
        c(proteinID_column, sequence_column, "Spectrum.File", quantification_column),
        c("ProteinName", "PeptideSequence", "Run", "Intensity"))
    input[["PeptideModifiedSequence"]] = paste(input[["PeptideSequence"]], 
                                               input[["Modifications"]], 
                                               sep = "_")
    input[, -which(colnames(input) == "PeptideSequence")]
}