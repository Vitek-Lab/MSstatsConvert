#' Clean raw Proteome Discoverer data
#' @inheritParams .cleanRawPDMSstats
.cleanRawPD = function(msstats_object, quantification_column, protein_id_column,
                       sequence_column, remove_shared) {
    if (getDataType(msstats_object) == "MSstatsTMT") {
        .cleanRawPDTMT(msstats_object, remove_shared, protein_id_column)
    } else {
        .cleanRawPDMSstats(msstats_object, quantification_column, 
                           protein_id_column, sequence_column, remove_shared)
    }
}

#' Clean raw PD output
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @param quantification_column chr, name of a column used for quantification.
#' @param protein_id_column chr, name of a column with protein IDs.
#' @param sequence_column chr, name of a column with peptide sequences.
#' @param remove_shared lgl, if TRUE, shared peptides will be removed.
#' @return data.table
#' @keywords internal
.cleanRawPDMSstats = function(msstats_object, quantification_column, protein_id_column,
                              sequence_column, remove_shared) {
    XProteins = NULL
    
    pd_input = getInputFile(msstats_object, "input")
    protein_id_column = .standardizeColnames(protein_id_column)
    sequence_column = .standardizeColnames(sequence_column)
    
    quantification_column = .findAvailable(c("Intensity", "Area"),
                                           colnames(pd_input),
                                           quantification_column)
    protein_id_column = .findAvailable(c("Protein.Accessions", 
                                         "MasterProteinAccessions",
                                         "ProteinGroupAccessions"),
                                       colnames(pd_input),
                                       protein_id_column)
    sequence_column = .findAvailable(c("Sequence", "AnnotatedSequence"), colnames(pd_input), 
                                     sequence_column)
    if (remove_shared) {
        pd_input = pd_input[XProteins == "1", ]
    }
    pd_cols = c(protein_id_column, sequence_column, 
                "Modifications", "Charge", "SpectrumFile", quantification_column)
    if (any(is.element(colnames(pd_input), "Fraction"))) {
        pd_cols = c(pd_cols, "Fraction")
    }
    pd_input = pd_input[, pd_cols, with = FALSE]
    colnames(pd_input) = .updateColnames(
        pd_input,
        c(protein_id_column, sequence_column, "SpectrumFile", quantification_column, "Charge"),
        c("ProteinName", "PeptideSequence", "Run", "Intensity", "PrecursorCharge"))
    pd_input[["PeptideSequence"]] = paste(pd_input[["PeptideSequence"]], 
                                          pd_input[["Modifications"]], 
                                          sep = "_")
    pd_input[, !(colnames(pd_input) == "Modifications"), with = FALSE]
}


#' Clean raw TMT data from Proteome Discoverer
#' @inheritParams .cleanRawPDMSstats
#' @importFrom data.table melt
#' @return `data.table`
#' @keywords internal
.cleanRawPDTMT = function(msstats_object, remove_shared = TRUE, protein_id_column = "ProteinAccessions") {
    pd_input = getInputFile(msstats_object, "input")
    protein_id_column = .standardizeColnames(protein_id_column)
    if (!is.element(protein_id_column, colnames(pd_input))) {
        protein_id_column = .findAvailable(c("ProteinAccessions", "MasterProteinAccessions"),
                                           colnames(pd_input), "ProteinAccessions")
    }
    if (protein_id_column == "ProteinAccessions") {
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
    
    pd_cols = intersect(c(protein_id_column, num_proteins, "AnnotatedSequence", "Charge", "IonsScore",
                          "SpectrumFile", "QuanInfo", "IsolationInterference", channels),
                        colnames(pd_input))
    pd_input = pd_input[, pd_cols, with = FALSE]
    colnames(pd_input) = .updateColnames(pd_input,
                                         c(protein_id_column, num_proteins, "AnnotatedSequence", "SpectrumFile"),
                                         c("ProteinName", "numProtein", "PeptideSequence", "Run"))
    pd_input$PSM = paste(pd_input$PeptideSequence, pd_input$Charge,
                         1:nrow(pd_input), sep = "_")
    pd_input = melt(pd_input, measure.vars = channels, 
                    id.vars = setdiff(colnames(pd_input), channels),
                    variable.name = "Channel", value.name = "Intensity")
    pd_input$Channel = .standardizeColnames(pd_input$Channel)
    pd_input = pd_input[(pd_input$ProteinName != "") & (!is.na(pd_input$ProteinName)), ]
    if (remove_shared) {
        if ("UNIQUE" %in% toupper(pd_input[["QuanInof"]])) {
            pd_input = pd_input[toupper(pd_input$QuanInfo) == 'UNIQUE', ]
        }
    }
    pd_input
}
