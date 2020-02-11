#' Import Spectronaut files
#' 
#' @param input name of Spectronaut output, which is long-format. ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity, F.ExcludedFromQuantification are required. Rows with F.ExcludedFromQuantification=True will be removed.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Spectronaut, use annotation=NULL (default). It will use the annotation information from input.
#' @param intensity 'PeakArea'(default) uses not normalized peak area. 'NormalizedPeakArea' uses peak area normalized by Spectronaut.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in EG.Qvalue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for EG.Qvalue. default is 0.01.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
SpectronauttoMSstatsFormat <- function(
  input, annotation = NULL, intensity = 'PeakArea', filter_with_Qvalue = TRUE,
  qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements="remove",
  removeProtein_with1Feature = FALSE, summaryforMultipleRows = max) {
  
  .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))
  .isLegalValue(intensity, legal_values = c("PeakArea", "NormalizedPeakArea"))

  .checkColumns("Input", 
                c("F.FrgLossType", "F.ExcludedFromQuantification",
                  "PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge",
                  "F.FrgIon", "R.FileName", "EG.Qvalue"), colnames(input))
  .checkColumns("Input", c("F.PeakArea", "F.NormalizedPeakArea"), 
                colnames(input), "optional")
  .checkColumns("Input", c("F.Charge", "F.FrgZ"), colnames(input), "optional")
  input = input[input[["F.FrgLossType"]] == "noloss", ]
  
  annotation = .makeAnnotation(
    annotation, 
    c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
    input, 
    c("R.FileName" = "Run", "R.Condition" = "Condition", "R.Replicate" = "BioReplicate")
  )  
  
  
  ## 2. use only 'F.ExcludedFromQuantification == False' : XIC quality
  if (is.logical(input$F.ExcludedFromQuantification)) {
    input$F.ExcludedFromQuantification <- factor(input$F.ExcludedFromQuantification)
    input$F.ExcludedFromQuantification <- factor(input$F.ExcludedFromQuantification,
                                                 labels = c('False', 'True'))
  }
  if (!all(unique(input$F.ExcludedFromQuantification) %in% c('False', 'True'))) {
    stop( paste("** Please check the column called F.ExcludedFromQuantification. Only False or True are allowed in this column."))
  }
  input <- input[input$F.ExcludedFromQuantification == 'False', ]
  
  f_charge_col = .findAvailable(c("F.Charge", "F.FrgZ"), colnames(input))
  pg_qval_col = findAvailable(c("PG.Qvalue"), colnames(input))
  input = .selectColumns(
    input, 
    c("PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge", "F.FrgIon", 
      f_charge_col, "R.FileName", "EG.Qvalue", pg_qval_col, paste0("F.", intensity)))
  colnames(input) = .updateColnames(
    input, 
    c("PG.ProteinGroups" = "ProteinName", "EG.ModifiedSequence" = "PeptideSequence",
      "FG.Charge" = "PrecursorCharge", "F.FrgIon" = "FragmentIon",
      f_charge_col = "ProductCharge", "R.FileName" = "Run", "EG.Qvalue" = "Qvalue",
      paste0("F.", intensity) = "Intensity"))

  ## 4. filter by Qvalue
  ## protein FDR
  if (is.element('PG.Qvalue', colnames(input))) {
    input[!is.na(input$PG.Qvalue) & input$PG.Qvalue > 0.01, "Intensity"] <- NA
    message('** Intensities with great than 0.01 in PG.Qvalue are replaced with NA.')
    input <- input[, -which(colnames(input) %in% 'PG.Qvalue')]
  }
  ## precursor qvalue
  if (filter_with_Qvalue) {
    if (!is.element(c('Qvalue'), colnames(input))) {
      stop('** EG.Qvalue column is needed in order to filter out by Qvalue. Please add EG.Qvalue column in the input.')
    } else {
      ## when qvalue > qvalue_cutoff, replace with zero for intensity
      input[!is.na(input$Qvalue) & input$Qvalue > qvalue_cutoff, "Intensity"] <- 0
      message(paste0('** Intensities with great than ', qvalue_cutoff, ' in EG.Qvalue are replaced with zero.'))
    }
  }
  
  ## 5. remove featuares with all na or zero
  ## some rows have all zero values across all MS runs. They should be removed.
  input$fea <- paste(input$PeptideSequence,
                     input$PrecursorCharge,
                     input$FragmentIon,
                     input$ProductCharge,
                     sep="_")
  inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
  count <- summarise(group_by(inputtmp, fea), length = length(Intensity))
  
  ## get feature with all NA or zeros
  getfea <- count[count$length > 0, 'fea']
  if (nrow(getfea) > 0) {
    nfea.remove <- length(unique(input$fea))-nrow(getfea)
    input <- input[which(input$fea %in% getfea$fea), ]
    message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
  }
  rm(inputtmp)
  
  input = .handleSharedPeptides(input, "ProteinName", "PeptideSequence",
                                remove_shared = useUniquePeptide)

  ##  7. remove features which has 1 or 2 measurements across runs
  if (fewMeasurements == "remove"){
    ## it is the same across experiments. # measurement per feature. 
    xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
    count_measure <- xtabs( ~fea, xtmp)
    remove_feature_name <- count_measure[count_measure < 3]
    if (length(remove_feature_name) > 0) {
      input <- input[-which(input$fea %in% names(remove_feature_name)), ]
      message(paste0('** ', length(remove_feature_name), ' features have 1 or 2 intensities across runs and are removed.'))
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
    input_w <- dcast( ProteinName + PeptideSequence + PrecursorCharge + FragmentIon + ProductCharge ~ Run, data=input, 
                      value.var='Intensity', 
                      fun.aggregate=summaryforMultipleRows, fill=NA_real_) 
    input <- melt(input_w, id=c('ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon', 'ProductCharge'))
    colnames(input)[which(colnames(input) %in% c('variable','value'))] <- c("Run","Intensity")
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
  } else {
    input <- input[, -which(colnames(input) %in% c('fea', 'Qvalue'))]
    message('** No multiple measurements in a feature and a run.')
  }
  
  ## 10. merge annotation
  input <- merge(input, annotation, all=TRUE)
  input <- data.frame("ProteinName" = input$ProteinName,
                      "PeptideSequence" = input$PeptideSequence,
                      "PrecursorCharge" = input$PrecursorCharge,
                      "FragmentIon" = input$FragmentIon,
                      "ProductCharge" = input$ProductCharge,
                      "IsotopeLabelType" = "L",
                      "Condition" = input$Condition,
                      "BioReplicate" = input$BioReplicate,
                      "Run" = input$Run,
                      "Intensity" = input$Intensity)
  
  .fixColumnTypes(input, factor_columns = "ProteinName")
}
