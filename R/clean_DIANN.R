#' Clean raw Diann files
#' @param msstats_object an object of class `MSstatsDIANNFiles`.
#' @param MBR True if analysis was done with match between runs
#' @return data.table
#' @keywords internal
.cleanRawDIANN = function(msstats_object, MBR = TRUE) {
  # read input and set data type
  dn_input = getInputFile(msstats_object, "input")
  dn_input = data.table::as.data.table(dn_input)
  
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
  
  # extract the intensities and fragment information
  dn_input = dn_input[, lapply(.SD, function(x) unlist(tstrsplit(x, ";"))),
                      .SDcols = c("FragmentQuantCorrected", "FragmentInfo"), 
                      by = setdiff(colnames(dn_input), c("FragmentInfo", "FragmentQuantCorrected"))]
  dn_input[, FragmentQuantCorrected := as.numeric(FragmentQuantCorrected)]
  dn_input[, FragmentIon := sub('\\^.*','',FragmentInfo)]
  dn_input[, ProductCharge := unlist(strsplit(FragmentInfo, split = '/'))[[1]], 
           by = FragmentInfo]
  dn_input[, ProductCharge := strtoi(sub('.*\\^','',ProductCharge))]
  # remove unnecessary fragment ions
  dn_input = dn_input[!grepl("NH3", FragmentIon), ]
  dn_input = dn_input[!grepl("H2O", FragmentIon), ]
  # remove rows with NA in the Fragment.Quant.Corrected and Fragment.Info columns
  dn_input = na.omit(dn_input, cols = c("FragmentInfo", "FragmentQuantCorrected"))
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
  protein.id = unique(dn_input$ProteinName)
  ## if the protein accession has multiple ids, there should be semicolon ;
  get.proid = protein.id[-grep(';', protein.id)]
  ## let's remove them
  dn_input = dn_input[ProteinName %in% get.proid, ]
  # log success and return results
  .logSuccess("DIANN", "clean")
  dn_input
}
