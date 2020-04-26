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

PDtoMSstatsFormat = function(
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
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide,
                                           feature_cols)
    input = merge(input, annotation, by = "Run")
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 "IsotopeLabelType" = "L"))
    input
}


#' Clean raw PD output
#' @param pd_input PD report or a path to it.
#' @param quantification_column chr, name of a column used for quantification.
#' @param proteinID_column chr, name of a column with protein IDs.
#' @param sequence_column chr, name of a column with peptide sequences.
#' @param filter_num_col lgl, if TRUE, shared peptides will be removed.
#' @param ... optional, additional parameters to data.table::fread.
#' @return data.table
#' @keywords internal
.cleanRawPD = function(pd_input, quantification_column, proteinID_column,
                       sequence_column, filter_num_col, ...) {
    pd_input = .getDataTable(pd_input, ...)
    colnames(pd_input) = .standardizeColnames(colnames(pd_input))
    quantification_column = .findAvailable(c("Intensity", "Area"),
                                           colnames(pd_input),
                                           quantification_column)
    # which.quantification = .isLegalValue(which.quantification, 
    #                                      legal_values = c("Intensity", "Area", "Precursor.Area",
    #                                                       "Precursor.Abundance"),
    #                                      message = "Please select a column to be used for quantified intensities among four options: ")
    
    proteinID_column = .findAvailable(c("Protein.Accessions", 
                                        "Master.Protein.Accessions"),
                                      colnames(pd_input),
                                      proteinID_column)
    # which.proteinid = .isLegalValue(which.proteinid, 
    #                                 legal_values = c("ProteinAccessions", 
    #                                                  "Master.Protein.Accessions",
    #                                                  "Protein.Group.Accessions"),
    #                                 message = "Please select a column to be used as protein IDs among three options: ")
    sequence_column = .findAvailable("Annotated.Sequence", colnames(pd_input), 
                                     sequence_column)
    # which.sequence = .isLegalValue(which.sequence, 
    #                                legal_values = c("Annotated.Sequence", "Sequence"),
    #                                message = "Please select peptide sequence column between two options: ")
    
    if(filter_num_col) {
        pd_input = pd_input[pd_input[["X..Proteins"]] == '1', ]
    }
    pd_cols = c(proteinID_column, "X..Proteins", sequence_column, 
                "Modifications", "Charge", "Spectrum.File", quantification_column)
    if (any(is.element(colnames(pd_input), "Fraction"))) {
        pd_cols = c(pd_cols, "Fraction")
    }
    pd_input = pd_input[, pd_cols, with = FALSE]
    colnames(pd_input) = .updateColnames(
        pd_input,
        c(proteinID_column, sequence_column, "Spectrum.File", quantification_column),
        c("ProteinName", "PeptideSequence", "Run", "Intensity"))
    pd_input[["PeptideModifiedSequence"]] = paste(pd_input[["PeptideSequence"]], 
                                                  pd_input[["Modifications"]], 
                                                  sep = "_")
    pd_input[, !(colnames(pd_input) == "PeptideSequence"), with = FALSE]
}


#' Convert Proteome Discoverer output to MSstatsTMT format.
#' 
#' @param input PD report or a path to it.
#' @param annotation annotation with Run, Fraction, TechRepMixture, Mixture, Channel, 
#' BioReplicate, Condition columns or a path to file. Refer to the example 'annotation' for the meaning of each column.
#' @param which.proteinid Use 'Proteins'(default) column for protein name. 'Leading.proteins' or 'Leading.razor.proteins' or 'Gene.names' can be used instead to get the protein ID with single protein. However, those can potentially have the shared peptides.
#' @param useUniquePeptide lgl, if TRUE (default) removes peptides that are assigned for more than one proteins. We assume to use unique peptide for each protein.
#' @param rmPSM_withMissing_withinRun lgl, if TRUE, will remove PSM with any missing value within each Run. Default is FALSE.
#' @param rmPSM_withfewMea_withinRun lgl, only for rmPSM_withMissing_withinRun = FALSE. TRUE (default) will remove the features that have 1 or 2 measurements within each Run.
#' @param rmProtein_with1Feature TRUE will remove the proteins which have only 1 peptide and charge. Defaut is FALSE.
#' @param summaryforMultipleRows sum(default) or max - when there are multiple measurements for certain feature in certain run, select the feature with the largest summation or maximal value.
#' @return data.table
#' @export
#' 
PDtoMSstatsTMTFormat <- function(
    input, annotation, which.proteinid = 'Protein.Accessions', 
    useNumProteinsColumn = TRUE, useUniquePeptide = TRUE, 
    rmPSM_withMissing_withinRun = FALSE, rmPSM_withfewMea_withinRun = TRUE, 
    rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum,
    use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
    .setMSstatsLogger(use_log_file, append, verbose)
    # annotation = .makeAnnotation()
    # TODO: checks!
    input = .cleanRawPDTMT(input, remove_shared = useUniquePeptide, protein_ID = which.proteinid)
    input = .filterExact(input, "numProtein", "1", TRUE, useNumProteinsColumn)
    feature_cols = c("PeptideSequence", "Charge")
    input = .removeMissingAllChannels(input, feature_cols)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .cleanByFeatureTMT(input, feature_cols,
                               summaryforMultipleRows, rmPSM_withfewMea_withinRun,
                               rmPSM_withMissing_withinRun)
    input = .mergeAnnotation(input, annotation)
    input = .handleSingleFeaturePerProtein(input, rmProtein_with1Feature, 
                                           feature_cols)
    input = .handleFractions(input, annotation)
    input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
              "TechRepMixture", "Run", "Channel", "Condition", "BioReplicate", "Intensity")] # unique?
}


#' Clean raw TMT data from Proteome Discoverer
#' @param pd_input PD report or a path to it.
#' @param remove_shared lgl, if TRUE, shared peptides will be removed.
#' @param protein_ID chr, names of a column with protein names.
#' @param ... optional, additional parameters for data.table::fread.
#' @importFrom data.table melt
#' @return data.table
#' @keywords internal
.cleanRawPDTMT = function(pd_input, remove_shared = TRUE, protein_ID = "ProteinAccessions", ...) {
    pd_input = .getDataTable(pd_input, ...)
    colnames(pd_input) = .standardizeColnames(colnames(pd_input))
    colnames(pd_input) = gsub(".", "", colnames(pd_input), fixed = TRUE)
    protein_ID = gsub(".", "", protein_ID, fixed = TRUE)
    if (!is.element(protein_ID, colnames(pd_input))) {
        protein_ID = .findAvailable(c("ProteinAccessions", "MasterProteinAccessions"),
                                    colnames(pd_input), "ProteinAccessions")
    }
    if (protein_ID == "ProteinAccessions") {
        num_proteins = "XProteins"
    } else {
        num_proteins = "XProteinsGroups"
    }
    
    channels = .getChannelColumns(colnames(pd_input), "Abundance")
    if (length(channels) == 0L) {
        msg = "There is no channel intensity column in the input data, which should start with 'Abundance'."
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    
    pd_cols = c(protein_ID, num_proteins, "AnnotatedSequence", "Charge", "IonsScore",
                "SpectrumFile", "QuanInfo", "IsolationInterference", channels)
    pd_input = pd_input[, pd_cols, with = FALSE]
    colnames(pd_input) = .updateColnames(pd_input,
                                         c(protein_ID, num_proteins, "AnnotatedSequence", "SpectrumFile"),
                                         c("ProteinName", "numProtein", "PeptideSequence", "Run"))
    pd_input$PSM = paste(pd_input$PeptideSequence, pd_input$Charge,
                         1:nrow(pd_input), sep = "_")
    pd_input = melt(pd_input, measure.vars = channels, 
                    id.vars = setdiff(colnames(pd_input), channels),
                    variable.name = "Channel", value.name = "Intensity")
    pd_input$Channel = gsub("Abundance", "", pd_input$Channel)
    pd_input$Channel = gsub("[:\\.]", "", pd_input$Channel)
    pd_input = pd_input[(pd_input$ProteinName != "") & (!is.na(pd_input$ProteinName)), ]
    if (remove_shared) {
        if ("UNIQUE" %in% toupper(pd_input[["QuanInof"]])) {
            pd_input = pd_input[toupper(pd_input$QuanInfo) == 'UNIQUE', ]
        }
    }
    pd_input
}
