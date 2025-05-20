#' Clean raw Metamorpheus files
#' @param msstats_object an object of class `MSstatsMetamorpheusFiles`.
#' @param MBR If TRUE, the function will include peaks detected by MBR
#' @param qvalue_cutoff The q-value cutoff for filtering peaks detected by MBR
#' @return data.table
#' @keywords internal
.cleanRawMetamorpheus = function(msstats_object, MBR = TRUE, qvalue_cutoff = 0.05) {
    metamorpheus_input = getInputFile(msstats_object, "input")
    metamorpheus_input = data.table::as.data.table(metamorpheus_input)
    # Remove peptide decoys 
    if ("DecoyPeptide" %in% colnames(metamorpheus_input)) {
        metamorpheus_input <- metamorpheus_input[!metamorpheus_input$DecoyPeptide, ]
    }
    # Remove peak decoys
    if ("RandomRT" %in% colnames(metamorpheus_input)) {
        metamorpheus_input <- metamorpheus_input[!metamorpheus_input$RandomRT, ]
    }
    if (MBR) {
        metamorpheus_input = metamorpheus_input[
            metamorpheus_input$PeakDetectionType == "MSMS" | 
                metamorpheus_input$`PIPQ-Value` <= qvalue_cutoff , ]
    } else {
        metamorpheus_input = metamorpheus_input[
            metamorpheus_input$PeakDetectionType == "MSMS", ]
    }
    req_cols = c('FileName', 'ProteinGroup', 'FullSequence', 
                 'PrecursorCharge', 'Peakintensity')
    metamorpheus_input = metamorpheus_input[, req_cols, with = FALSE]
    data.table::setnames(
        metamorpheus_input, 
        c("ProteinGroup", "FullSequence", "PrecursorCharge", "FileName", "Peakintensity"),
        c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run", "Intensity"), 
        skip_absent = TRUE)
    .logSuccess("Metamorpheus", "clean")
    metamorpheus_input
}