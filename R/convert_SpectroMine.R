#' Clean raw SpectroMine TMT data
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @importFrom data.table melt
#' @return `data.table`
#' @keywords internal
.cleanRawSpectroMineTMT = function(msstats_object) {
  PSM = PeptideSequence = PrecursorCharge = ProteinName = NULL
  
  sm_input = getInputFile(msstats_object, "input")
  channels = .getChannelColumns(colnames(sm_input), "PSM", "Raw")
  if (length(channels) == 0L) {
    msg = paste("There is no channel intensity column in the input data,", 
                "which should start with 'PSM' and end with 'Raw'.")
    getOption("MSstatsMsg")("ERROR", msg)
    stop(msg)
  }
  sm_input = sm_input[, c("PGProteinAccessions", "PMoleculeID", "PPCharge",
                          "PGQValue", "PSMQvalue", "RFileName", channels),
                      with = FALSE]
  data.table::setnames(sm_input, 
                       c("PGProteinAccessions", "PMoleculeID", 
                         "PPCharge", "PSMQvalue", "RFileName"),
                       c("ProteinName", "PeptideSequence", "PrecursorCharge",
                         "Qvalue", "Run"))
  sm_input = sm_input[(ProteinName != "") & (!is.na(ProteinName)), ]
  sm_input[, PSM := paste(PeptideSequence, PrecursorCharge, 
                          1:nrow(sm_input), sep = "_")]
  sm_input = melt(sm_input, measure.vars = channels,
                  id.vars = setdiff(colnames(sm_input), channels),
                  variable.name = "Channel", value.name = "Intensity")
  sm_input$Channel = gsub("PSM", "", sm_input$Channel)
  sm_input$Channel = gsub("Raw", "", sm_input$Channel)
  sm_input$Channel = gsub(".", "", sm_input$Channel, fixed = TRUE)
  sm_input$Intensity = ifelse(sm_input$Intensity == 0, NA, sm_input$Intensity)
  sm_input
}
