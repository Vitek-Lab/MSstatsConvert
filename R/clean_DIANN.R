#' Clean raw Diann files
#' @param msstats_object an object of class `MSstatsDIANNFiles`.
#' @param MBR True if analysis was done with match between runs
#' @return data.table
#' @importFrom stats na.omit
#' @keywords internal
.cleanRawDIANN = function(msstats_object, MBR = TRUE) {
  dn_input = getInputFile(msstats_object, "input")
  dn_input = data.table::as.data.table(dn_input)
  
  if (!is.element("PrecursorMz", colnames(dn_input))) {
    dn_input[, PrecursorMz := NA]
  }
  if (!is.element('FragmentInfo', colnames(dn_input))) {
    dn_input[, FragmentInfo := NA]
  }
  req_cols = c('ProteinNames', 'StrippedSequence', 
               'ModifiedSequence', 'PrecursorCharge',
               'FragmentQuantCorrected', 'QValue', 
               'PrecursorMz', 'FragmentInfo', 'Run')
  if (MBR) {
    req_cols = c(req_cols, c('LibQValue', 'LibPGQValue'))
  } else{
    req_cols = c(req_cols, c('GlobalQValue', 'GlobalPGQValue'))
  }
  dn_input = dn_input[, req_cols, with = FALSE]
  dn_input = dn_input[, lapply(.SD, function(x) unlist(tstrsplit(x, ";"))),
                      .SDcols = c("FragmentQuantCorrected", "FragmentInfo"), 
                      by = setdiff(colnames(dn_input), c("FragmentInfo", "FragmentQuantCorrected"))]
  if (all(is.na(dn_input[["FragmentInfo"]]))) {
    dn_input[, FragmentInfo := paste0("Frag", 1:.N),
             by = c("ProteinNames", "ModifiedSequence", "PrecursorCharge", "Run")]
  }
  dn_input[, FragmentQuantCorrected := as.numeric(FragmentQuantCorrected)]
  dn_input[, FragmentIon := sub('\\^\\.\\*', '', FragmentInfo)]
  if (any(grepl("/", dn_input$FragmentInfo))) {
    dn_input[, ProductCharge := unlist(strsplit(FragmentInfo, split = "/"))[[1]], by = FragmentInfo]
    dn_input[, ProductCharge := strtoi(sub("\\.\\*\\^", "", ProductCharge))]
  } else {
    dn_input[, ProductCharge := 1]
  }
  dn_input = dn_input[!grepl("NH3", FragmentIon), ]
  dn_input = dn_input[!grepl("H2O", FragmentIon), ]
  dn_input = na.omit(dn_input, cols = "FragmentQuantCorrected")
  data.table::setnames(dn_input, old = c('ProteinNames', 'StrippedSequence', 
                                         'ModifiedSequence','PrecursorCharge',
                                         'FragmentQuantCorrected', 'QValue', 
                                         'PrecursorMz', 'FragmentIon','Run', 
                                         'ProductCharge'),
                       new = c('ProteinName', 'PeptideSequence', 
                               'PeptideModifiedSequence','PrecursorCharge',
                               'Intensity', 'DetectionQValue', 
                               'PrecursorMz', 'FragmentIon','Run',
                               'ProductCharge'),
                       skip_absent = TRUE)
  dn_input[, PeptideSequence := NULL]
  setnames(dn_input, "PeptideModifiedSequence", "PeptideSequence")
  .logSuccess("DIANN", "clean")
  dn_input
}
