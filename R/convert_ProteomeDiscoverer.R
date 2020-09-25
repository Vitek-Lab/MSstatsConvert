#' Clean raw Proteome Discoverer data
#' @inheritParams .cleanRawPDTMT
#' @inheritParams .cleanRawPDMSstats
.cleanRawPD = function(msstats_object, quantification_column, protein_id_column,
                       sequence_column, remove_shared, remove_protein_groups = TRUE,
                       intensity_columns_regexp = "Abundance") {
    if (getDataType(msstats_object) == "MSstatsTMT") {
        .cleanRawPDTMT(msstats_object, remove_shared, remove_protein_groups, 
                       protein_id_column, intensity_columns_regexp)
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
    quantification_column = .standardizeColnames(quantification_column)
    
    quantification_column = .findAvailable(c("Intensity", "Area", "PrecursorArea"),
                                           colnames(pd_input),
                                           quantification_column)
    protein_id_column = .findAvailable(c("ProteinAccessions", 
                                         "MasterProteinAccessions",
                                         "ProteinGroupAccessions"),
                                       colnames(pd_input),
                                       protein_id_column)
    sequence_column = .findAvailable(c("Sequence", "AnnotatedSequence"), 
                                     colnames(pd_input), 
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
    data.table::setnames(
        pd_input,
        c(protein_id_column, sequence_column, "SpectrumFile", quantification_column, "Charge"),
        c("ProteinName", "PeptideSequence", "Run", "Intensity", "PrecursorCharge"),
        skip_absent = TRUE)
    pd_input[["PeptideSequence"]] = paste(pd_input[["PeptideSequence"]], 
                                          pd_input[["Modifications"]], 
                                          sep = "_")
    pd_input[, !(colnames(pd_input) == "Modifications"), with = FALSE]
}


#' Clean raw TMT data from Proteome Discoverer
#' @inheritParams .cleanRawPDMSstats
#' @param remove_protein_groups if TRUE, proteins with numProteins > 1 will be removed.
#' @param intensity_columns_regexp regular expressions that defines intensity columns.
#' Defaults to "Abundance", which means that columns that contain the word "Abundance"
#' will be treated as corresponding to intensities for different channels.
#' @importFrom data.table melt
#' @return `data.table`
#' @keywords internal
.cleanRawPDTMT = function(msstats_object, remove_shared = TRUE, 
                          remove_protein_groups = TRUE,
                          protein_id_column = "ProteinAccessions",
                          intensity_columns_regexp = "Abundance") {
    ProteinName = numProtein = QuanInfo = NULL
    
    pd_input = getInputFile(msstats_object, "input")
    protein_id_column = .standardizeColnames(protein_id_column)
    if (!is.element(protein_id_column, colnames(pd_input))) {
        protein_id_column = .findAvailable(c("ProteinAccessions", "MasterProteinAccessions"),
                                           colnames(pd_input), "ProteinAccessions")
    }
    if (protein_id_column == "ProteinAccessions") {
        num_proteins = "XProteins"
    } else {
        num_proteins = "XProteinGroups"
    }
    
    channels = .getChannelColumns(colnames(pd_input), intensity_columns_regexp)
    if (length(channels) == 0L) {
        msg = paste("There is no channel intensity column in the input data,",
                    "which should start with the string provided in the",
                    "intensity_columns_regexp parameter (default: 'Abundance')")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    
    pd_cols = intersect(c(protein_id_column, num_proteins, "AnnotatedSequence", "Charge", "PrecursorCharge",
                          "IonsScore", "SpectrumFile", "QuanInfo", "IsolationInterference", channels),
                        colnames(pd_input))
    pd_input = pd_input[, pd_cols, with = FALSE]
    data.table::setnames(pd_input,
                         c(protein_id_column, num_proteins, "AnnotatedSequence", "SpectrumFile", "Charge"),
                         c("ProteinName", "numProtein", "PeptideSequence", "Run", "PrecursorCharge"),
                         skip_absent = TRUE)
    pd_input = unique(pd_input)
    pd_input$PSM = paste(pd_input$PeptideSequence, pd_input$PrecursorCharge,
                         1:nrow(pd_input), sep = "_")
    pd_input = melt(pd_input, measure.vars = channels, 
                    id.vars = setdiff(colnames(pd_input), channels),
                    variable.name = "Channel", value.name = "Intensity")
    pd_input$Channel = .standardizeColnames(pd_input$Channel)
    pd_input$Channel = gsub(intensity_columns_regexp, "", pd_input$Channel)
    pd_input$Channel = gsub(":", "", pd_input$Channel)
    pd_input$Intensity = ifelse(pd_input$Intensity == 0, NA, pd_input$Intensity)
    pd_input = pd_input[(ProteinName != "") & (!is.na(ProteinName)), ]
    if (remove_protein_groups) {
        pd_input = pd_input[numProtein == 1]
    }
    if (remove_shared) {
        if ("UNIQUE" %in% toupper(pd_input[["QuanInfo"]])) {
            pd_input = pd_input[toupper(QuanInfo) == 'UNIQUE', ]
        }
    }
    pd_input
}
