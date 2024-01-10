#' Clean raw Metamorpheus files
#' @param msstats_object an object of class `MSstatsMetamorpheusFiles`.
#' @return data.table
#' @keywords internal
.cleanRawMetamorpheus = function(msstats_object) {
    metamorpheus_input = getInputFile(msstats_object, "input")
    metamorpheus_input = data.table::as.data.table(metamorpheus_input)
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