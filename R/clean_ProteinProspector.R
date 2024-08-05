#' Clean raw Protein Prospector data
#' @param msstats_object object that inherits from MSstatsInputFiles class.
#' @return data.table
#' @keywords internal
#' @noRd
.cleanRawProteinProspector = function(msstats_object) {
    protein_prospector_input = getInputFile(msstats_object, "input")
    protein_prospector_input = 
        data.table::as.data.table(protein_prospector_input)
    channels = .getChannelColumns(
        colnames(protein_prospector_input), "Int")
    req_cols = c('AccX', 'z', 'DBPeptide', 'Fraction', channels)
    protein_prospector_input = protein_prospector_input[, req_cols, with = FALSE]
    data.table::setnames(
        protein_prospector_input, 
        c("AccX", "DBPeptide", "z", "Fraction"),
        c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run"), 
        skip_absent = TRUE)
    protein_prospector_input$PSM = paste(
        protein_prospector_input$PeptideSequence, 
        protein_prospector_input$PrecursorCharge,
        1:nrow(protein_prospector_input), sep = "_"
    )
    
    protein_prospector_input = melt(protein_prospector_input, 
                                    measure.vars = channels, 
                                    id.vars = setdiff(
                                        colnames(protein_prospector_input), 
                                        channels
                                    ),
                                    variable.name = "Channel", 
                                    value.name = "Intensity"
                                )
    protein_prospector_input$Channel = .standardizeColnames(
        protein_prospector_input$Channel
    )
    
    .logSuccess("ProteinProspector", "clean")
    protein_prospector_input
}