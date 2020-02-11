#' Import Proteome Discoverer files
#' 
#' @inheritParams .documentFunction
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, 
#' Run information. 'Run' will be matched with 'Spectrum.File'.
#' @param useNumProteinsColumn TRUE removes peptides which have more than 1 in # Proteins column of PD output.
#' @param which.quantification Use 'Precursor.Area'(default) column for quantified intensities. 'Intensity' or 'Area' can be used instead.
#' @param which.proteinid Use 'Protein.Accessions'(default) column for protein name. 'Master.Protein.Accessions' can be used instead.
#' @param which.sequence Use 'Sequence'(default) column for peptide sequence. 'Annotated.Sequence' can be used instead.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
PDtoMSstatsFormat <- function(input,
                              annotation,
                              useNumProteinsColumn=FALSE,
                              useUniquePeptide=TRUE,
                              summaryforMultipleRows=max,
                              fewMeasurements="remove",
                              removeOxidationMpeptides=FALSE,
                              removeProtein_with1Peptide=FALSE,
                              which.quantification = 'Precursor.Area',
                              which.proteinid = 'Protein.Group.Accessions',
                              which.sequence = 'Sequence'){
    
    .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))

    annotation = .makeAnnotation(annotation, 
                                 c("Run" = "Run", "Condition" = "Condition",
                                   "BioReplicate" = "BioReplicate"))    
    ## 0.1. which intensity : Precursor.Area vs. Intensity vs Area
    ## 2017.01.11 : use 'Precursor.Area' instead of 'Intensity'
    ## default : Precursor.Area
    .isLegalValue(which.quantification, 
                  legal_values = c("Intensity", "Area", "Precursor.Area",
                                   "Precursor.Abundance"),
                  message = "Please select a column to be used for quantified intensities among four options: ")
    which.quant = which.quantification
    
    if (which.quant == 'Intensity' & !is.element('Intensity', colnames(input))) {
        ## then that is because, input came from different version
        which.quant <- 'Precursor.Area'
        message('** Use Precursor.Area instead of Intensity.')
    }
    
    if (which.quant == 'Area' & !is.element('Area', colnames(input))) {
        ## then that is because, input come from different version
        which.quant <- 'Precursor.Area'
        message('** Use Precursor.Area instead of Intensity.')
    }
    if (which.quant == 'Precursor.Area' & !is.element('Precursor.Area', colnames(input))) {
        ## then that is because, input come from different version 
        stop('** Please select which columns should be used for quantified intensities, among two options (Intensity or Area).')
    }
    if (!is.element(which.quant, colnames(input))) {
        stop('** Please select which columns should be used for quantified intensities, among three options (Intensity, Area, Precursor.Area).')
    }
    
    ## 0.2. which protein id : Protein Accessions vs Master Protein Accesisions
    ## default : Protein Accessions
    .isLegalValue(which.proteinid, 
                  legal_values = c("ProteinAccessions", 
                                   "Master.Protein.Accessions",
                                   "Protein.Group.Accessions"),
                  message = "Please select a column to be used as protein IDs among three options: ")
    which.pro = which.proteinid
    # FIXME: remove unnecessary renaming
    if (which.pro == 'Protein.Accessions' & !is.element('Protein.Accessions', colnames(input))) {
        which.pro <- 'Protein.Group.Accessions'
        message('** Use Protein.Group.Accessions instead of Protein.Accessions.')
    }
    if (which.pro == 'Master.Protein.Accessions' & !is.element('Master.Protein.Accessions', colnames(input))) {
        which.pro <- 'Protein.Group.Accessions'
        message('** Use Protein.Group.Accessions instead of Master.Protein.Accessions.')
    }
    if (which.pro == 'Protein.Group.Accessions' & !is.element('Protein.Group.Accessions', colnames(input))) {
        ## then that is because, input come from different version 
        stop('** Please select which columns should be used for protein ids, among two options (Protein.Accessions or Master.Protein.Accessions).')
    }
    if (!is.element(which.pro, colnames(input))) {
        stop('** Please select which columns should be used for protein ids, among three options (Protein.Accessions, Master.Protein.Accessions, Protein.Group.Accessions).')
    }
    
    ################################################
    ## 0.3. which sequence : Sequence vs Annotated.Sequence
    ################################################
    ## default : Sequence
    .isLegalValue(which.sequence, 
                  legal_values = c("Annotated.Sequence", "Sequence"),
                  message = "Please select peptide sequence column between two options: ")
    which.seq = which.sequence
    # FIXME: remove unncessary re-naming
    
    if (which.seq == 'Annotated.Sequence' & !is.element('Annotated.Sequence', colnames(input))) {
        which.seq <- 'Sequence'
        message('** Use Sequence instead of Annotated.Sequence.')
    }
    if (!is.element(which.seq, colnames(input))) {
        stop('** Please select which columns should be used for peptide sequence, between twp options (Sequence or Annotated.Sequence).')
    }
    
    ## 1. get subset of columns
    input = input[, c(which.pro, "X..Proteins", which.seq, "Modifications", 
                      "Charge", "Spectrum.File", which.quant)]
    colnames(input) = .updateColnames(
        input, c(which.pro = "ProteinName", "X..Proteins" = "numProtein",
                 which.seq = "PeptideSequence", "Spectrum.File" = "Run", 
                 which.quant = "Intensity"))
    ## 2. remove peptides which are used in more than one protein
    ## we assume to use unique peptide
    if (useNumProteinsColumn) {
        ## remove rows with #proteins is not 1
        input <- input[input$numProtein == '1', ]
        message('** Rows with #Proteins, which are not equal to 1, are removed.')
    }
    
    input = .handleSharedPeptides(input, "ProteinName", "PeptideSequence",
                                  remove_shared = useUniquePeptide)
    ## 3. remove the peptides including oxidation (M) sequence
    if (removeOxidationMpeptides) {
        remove_m_sequence <- unique(input[grep("Oxidation", input$Modifications), "Modifications"])
        if (length(remove_m_sequence) > 0) {
            input <- input[-which(input$Modifications %in% remove_m_sequence), ]
        }
        message('Peptides including oxidation(M) in the Modifications are removed.')
    }
    
    ## 4. remove multiple measurements per feature and run
    ## maximum or sum up abundances among intensities for identical features within one run
    input_sub <- dcast( ProteinName + PeptideSequence + Modifications + Charge ~ Run, data=input, 
                        value.var='Intensity', 
                        fun.aggregate=summaryforMultipleRows, 
                        fill=NA_real_) 
    ## reformat for long format
    input_sub <- melt(input_sub, id=c('ProteinName', 'PeptideSequence', 'Modifications', 'Charge'))
    colnames(input_sub)[which(colnames(input_sub) %in% c('variable','value'))] <- c("Run", "Intensity")
    message('** Multiple measurements in a feature and a run are summarized by summaryforMultipleRows.')
    input <- input_sub
    ## 5. add annotation
    noruninfo <- setdiff(unique(input$Run), unique(annotation$Run))
    if (length(noruninfo) > 0) {
        stop(paste('** Annotation for Run :', 
                   paste(noruninfo, collapse = ', '), 
                   ' are needed. Please update them in annotation file.') )
    }
    input <- merge(input, annotation, by="Run", all=TRUE)
    ## add other required information
    input$FragmentIon <- NA
    input$ProductCharge <- NA
    input$IsotopeLabelType <- "L"
    input$PeptideModifiedSequence <- paste(input$PeptideSequence, input$Modifications, sep="_")
    input.final <- data.frame("ProteinName" = input$ProteinName,
                              "PeptideModifiedSequence" = input$PeptideModifiedSequence,
                              "PrecursorCharge" = input$Charge,
                              "FragmentIon" = input$FragmentIon,
                              "ProductCharge" = input$ProductCharge,
                              "IsotopeLabelType" = input$IsotopeLabelType,
                              "Condition" = input$Condition,
                              "BioReplicate" = input$BioReplicate,
                              "Run" = input$Run,
                              "Intensity" = input$Intensity)
    if (any(is.element(colnames(input), 'Fraction'))) {
        input.final <- data.frame(input.final,
                                  "Fraction" = input$Fraction)
    }
    input <- input.final
    rm(input.final)
    
    ## 6. remove featuares with all na or zero
    ## some rows have all zero values across all MS runs. They should be removed.
    input$fea <- paste(input$PeptideModifiedSequence,
                       input$PrecursorCharge,
                       input$FragmentIon,
                       input$ProductCharge,
                       sep="_")
    inputtmp <- input[!is.na(input$Intensity) & input$Intensity > 1, ]
    count <- inputtmp %>% group_by(fea) %>% summarise(length=length(Intensity))
    ## get feature with all NA or zeros
    getfea <- count[count$length > 0, 'fea']
    if (nrow(getfea) > 0) {
        nfea.remove <- length(unique(input$fea))-nrow(getfea)
        input <- input[which(input$fea %in% getfea$fea), ]
        message(paste0('** ', nfea.remove, ' features have all NAs or zero intensity values and are removed.'))
    }
    rm(inputtmp)
    input <- input[, -which(colnames(input) %in% c('fea'))]
    
    ##  7. remove features which has 1 or 2 measurements across runs
    if (fewMeasurements=="remove") {
        ## it is the same across experiments. # measurement per feature. 
        xtmp <- input[!is.na(input$Intensity) & input$Intensity > 0, ]
        xtmp$feature <- paste(xtmp$PeptideModifiedSequence, xtmp$PrecursorCharge, sep="_")
        count_measure <- xtabs( ~feature, xtmp)
        remove_feature_name <- count_measure[count_measure < 3]
        input$feature <- paste(input$PeptideModifiedSequence, input$PrecursorCharge, sep="_")
        if (length(remove_feature_name) > 0) {
            input <- input[-which(input$feature %in% names(remove_feature_name)), ]
        }
        input <- input[, -which(colnames(input) %in% c('feature'))]
    }
    
    ##  8. remove proteins with only one peptide and charge per protein
    if (removeProtein_with1Peptide) {
        ## remove protein which has only one peptide
        input$feature <- paste(input$PeptideModifiedSequence, 
                               input$PrecursorCharge, 
                               input$FragmentIon, 
                               input$ProductCharge, 
                               sep="_")
        tmp <- unique(input[, c("ProteinName", 'feature')])
        tmp$ProteinName <- factor(tmp$ProteinName)
        count <- xtabs( ~ ProteinName, data=tmp)
        lengthtotalprotein <- length(count)
        removepro <- names(count[count <= 1])
        if (length(removepro) > 0) {
            input <- input[-which(input$ProteinName %in% removepro), ]
            message(paste0("** ", length(removepro), 
                           ' proteins, which have only one feature in a protein, are removed among ', 
                           lengthtotalprotein, ' proteins.'))
        }
        input <- input[, -which(colnames(input) %in% c('feature'))]
    }
    
    return(input)
}
