#' Clean raw OpenSWATH files
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @return data.table
#' @keywords internal
.cleanRawOpenSWATH = function(msstats_object) {
  PeptideSequence = FragmentIon = Intensity = NULL
  
  os_input = getInputFile(msstats_object, "input")
  os_input = os_input[, c("ProteinName", "FullPeptideName", "Charge", 
                          "filename", "aggr_Fragment_Annotation", "aggr_Peak_Area",
                          "m_score", "decoy"), with = FALSE]
  data.table::setnames(
    os_input,
    c("FullPeptideName", "Charge", "filename", "aggr_Fragment_Annotation", "aggr_Peak_Area"),
    c("PeptideSequence", "PrecursorCharge", "Run", "FragmentIon", "Intensity"),
    skip_absent = TRUE)
  os_input = os_input[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                          "Run", "FragmentIon", "Intensity",
                          "m_score", "decoy"), with = FALSE]
  os_input$Intensity = as.character(os_input$Intensity)
  os_input = os_input[, lapply(.SD, 
                               function(x) unlist(tstrsplit(x, split = ";", fixed = TRUE))),
                      by = c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run", "m_score", "decoy"),
                      .SDcols = c("FragmentIon", "Intensity")]
  os_input[, c("PeptideSequence", "FragmentIon")] = os_input[, lapply(list(PeptideSequence, FragmentIon),
                                                                   function(x) gsub(":", "_", x))]
  os_input[["Intensity"]] = as.numeric(as.character(os_input[["Intensity"]]))
  os_input$Intensity[os_input$Intensity < 1] = NA
  os_input
}
