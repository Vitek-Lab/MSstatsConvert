#' Import OpenSWATH files
#' 
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from OpenSWATH, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. 
#' Run should be the same as filename.
#' @param filter_with_mscore TRUE(default) will filter out the features that have greater than mscore_cutoff in m_score column. Those features will be removed.
#' @param mscore_cutoff Cutoff for m_score. default is 0.01.
#'  
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 

OpenSWATHtoMSstatsFormat <- function(
  input, annotation = NULL, filter_with_mscore = TRUE, mscore_cutoff = 0.01,
  useUniquePeptide = TRUE, fewMeasurements = "remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max) {
  
  .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))
  
  input = .selectColumns(
    input, 
    c("ProteinName", "FullPeptideName", "Charge", "Sequence", "decoy", "m_score",
      "aggr_Fragment_Annotation", "aggr_Peak_Area", "filename")
  )
  # TODO: only choose m_score if filter_with_mscore = TRUE
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate")
  )
  
  ## 2. remove the decoys
  if (length(unique(input$decoy)) == 2) {
    ## decoy=1 means the decoys.
    input <- input[input$decoy == 0, ]
    input <- input[, -which(colnames(input) %in% 'decoy')]
    message("** Remove the decoys.")
  }
  
  input = .handleFiltering(input, "m_score", mscore_cutoff, "smaller", "remove", 
                           NULL, TRUE, filter_with_mscore)

  ## 4. Make required long format - disaggregate : one row for each transition
  ## The columns "aggr_Fragment_Annotation" : separate by ';' and "aggr_Peak_Area" : separate by ';' 
  ## are disaggregated into the new columns "FragmentIon" and "Intensity". 
  input <- separate_rows(input, aggr_Fragment_Annotation, aggr_Peak_Area, sep = "[;]")
  
  # TODO: add difference between peptide name and full peptide name to documentation
  colnames(input) = .updateColnames(
    input, c("aggr_Fragment_Annotation" = "FragmentIon",
             "aggr_Peak_Area" = "Intensity", "FullPeptideName" = "PeptideSequence",
             "Charge" = "PrecursorCharge", "filename" = "Run"))
  input = .removeColumns(input, "Sequence")

  ## Unimod Identifier should be replaced from ":" to "_".
  input$PeptideSequence <- gsub(':', '_', input$PeptideSequence)
  input$FragmentIon <- gsub(':', '_', input$FragmentIon)
  
  ## there are incomplete rows, which potentially NA
  ## if negative and 0 values should be replaced with NA
  input[input$Intensity < 1, "Intensity"] <- 0
  
  ## 5. remove featuares with all na or zero
  ## some rows have all zero values across all MS runs. They should be removed.
  input$fea <- paste(input$PeptideSequence,
                     input$PrecursorCharge,
                     input$FragmentIon,
                     #input$ProductCharge,
                     sep="_")
  inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
  count <- summarise(group_by(inputtmp, fea), length=length(Intensity))
  ## get feature with all NA or zeros
  getfea <- count[count$length > 0, 'fea']
  if (nrow(getfea) > 0) {
    nfea.remove <- length(unique(input$fea)) - nrow(getfea)
    input <- input[which(input$fea %in% getfea$fea), ]
    message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
  }
  rm(inputtmp)
  
  input = .handleSharedPeptides(input, "ProteinName", "PeptideSequence",
                                remove_shared = useUniquePeptide)
  ##  7. remove features which has 1 or 2 measurements across runs
  if (fewMeasurements=="remove") {
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
  
  ## 8. remove proteins with only one peptide and charge per protein
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
  
  ## 9. remove multiple measurements per feature and run
  count <- aggregate(Intensity ~ fea, data=input, FUN=length)
  ## if any feature has more number of total MS runs, 
  if (any(unique(count$Intensity) > length(unique(input$Run)))) {
    ## maximum or sum up abundances among intensities for identical features within one run
    input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon ~ Run, data=input, 
                      value.var='Intensity', 
                      fun.aggregate=summaryforMultipleRows, fill='0') 
    ## reformat for long format
    input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  } else {
    ## still need to fill incomplete rows
    input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon ~ Run, data=input, 
                      value.var='Intensity', 
                      fill='0') 
    ## reformat for long format
    input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    message('** No multiple measurements in a feature and a run.')
  }
  
  ## 10. class of intensity is character, change it as numeric
  input$Intensity <- as.numeric(input$Intensity)
  
  ## 11. merge annotation
  input <- merge(input, annotinfo, by='Run', all=TRUE)
  
  ## fill in extra columns
  input.final <- data.frame("ProteinName" = input$ProteinName,
                            "PeptideSequence" = input$PeptideSequence,
                            "PrecursorCharge" = input$PrecursorCharge,
                            "FragmentIon" = input$FragmentIon,
                            "ProductCharge" = NA,
                            "IsotopeLabelType" = "L",
                            "Condition" = input$Condition,
                            "BioReplicate" = input$BioReplicate,
                            "Run" = input$Run,
                            "Intensity" = input$Intensity)
  input <- input.final
  input$ProteinName <- factor(input$ProteinName)
  rm(input.final)
  return(input)
}
