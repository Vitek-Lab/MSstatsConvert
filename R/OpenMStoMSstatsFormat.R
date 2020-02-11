#' Import OpenMS files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenMS, which includes feature(peptide ion)-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 

OpenMStoMSstatsFormat <- function(
  input, annotation = NULL, useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max) {
  
  .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))
  input = .selectColumns(input, 
                         c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                           "FragmentIon", "ProductCharge", "IsotopeLabelType",
                           "Condition", "BioReplicate", "Run", "Intensity"))
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
    input
  )

  ## 2. remove features with all na or zero
  ## some rows have all zero values across all MS runs. They should be removed.
  input$fea <- paste(input$PeptideSequence,
                     input$PrecursorCharge,
                     input$FragmentIon,
                     input$ProductCharge,
                     sep="_")
  inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
  count <- summarise(group_by(inputtmp, fea), length=length(Intensity))
  ## get feature with all NA or zeros
  getfea <- count[count$length > 0, 'fea']
  if (nrow(getfea) > 0) {
    nfea.remove <- length(unique(input$fea)) - nrow(getfea)
    input <- input[which(input$fea %in% getfea$fea), ]
    message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
  } else {
    stop(message('No intensity is available. Please check the input.'))
  }
  rm(inputtmp)
  
  input = .handleSharedPeptides(input, "ProteinName", "PeptideSequence",
                                remove_shared = useUniquePeptide)
  
  ##  4. remove features which has 1 or 2 measurements across runs
  if (fewMeasurements == "remove") {
    ## it is the same across experiments. # measurement per feature. 
    xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
    count_measure <- xtabs( ~fea, xtmp)
    remove_feature_name <- count_measure[count_measure < 3]
    if (length(remove_feature_name) > 0) {
      input <- input[-which(input$fea %in% names(remove_feature_name)), ]
      message(paste0('** ', length(remove_feature_name), 
                     ' features have 1 or 2 intensities across runs and are removed.'))
    }
  }
  
  ## 5. remove proteins with only one peptide and charge per protein
  if (removeProtein_with1Feature) {
    ## remove protein which has only one peptide
    tmp <- unique(input[, c("ProteinName", 'fea')])
    tmp$ProteinName <- factor(tmp$ProteinName)
    count <- xtabs( ~ ProteinName, data=tmp)
    lengthtotalprotein <- length(count)
    removepro <- names(count[count <= 1])
    if (length(removepro) > 0) {
      input <- input[-which(input$ProteinName %in% removepro), ]
      message(paste0("** ", length(removepro), 
                     ' proteins, which have only one feature in a protein, are removed among ', 
                     lengthtotalprotein, ' proteins.'))
    } else {
      message("** All proteins have at least two features.")
    }
  }
  
  ## 6. remove multiple measurements per feature and run
  count <- aggregate(Intensity ~ fea, data=input, FUN=length)
  ## if any feature has more number of total MS runs, 
  if (any(unique(count$Intensity) > length(unique(input$Run)))) {
    ## maximum or sum up abundances among intensities for identical features within one run
    input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon ~ Run, data=input, 
                      value.var='Intensity', 
                      fun.aggregate=summaryforMultipleRows, fill='NA') 
    ## reformat for long format
    input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  } else {
    ## still need to fill incomplete rows
    input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                      value.var='Intensity', 
                      fill='NA') 
    ## reformat for long format
    input <- melt(input_w, 
                  id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon', 'ProductCharge'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    message('** No multiple measurements in a feature and a run.')
  }
  
  input = merge(input, annotation, by = "Run", all = TRUE)
  input[["IsotopeLabelType"]] = "L"
  # TODO: check if the previous code did the intended thing
  .fixColumnTypes(input, factor_columns = "ProteinName", 
                  numeric_columns = "Intensity")
}
