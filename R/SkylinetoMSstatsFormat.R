#' Import Skyline files
#'
#' @inheritParams .documentFunction
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek
#' 
#' @export
#' 
SkylinetoMSstatsFormat <- function(
    input, annotation = NULL, removeiRT = TRUE, filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, useUniquePeptide = TRUE, fewMeasurements = "remove",
    removeOxidationMpeptides = FALSE, removeProtein_with1Feature = FALSE) {
    
    fewMeasurements = .isLegalValue(fewMeasurements, legal_values = c("remove", "keep"))
    
    input = .cleanRawSkyline(input)
    annotation = .makeAnnotation(
        annotation, 
        c("Run" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate"),
        input
    )
    input = .handleDecoyProteins(input, "ProteinName", c("DECOY", "Decoys"), FALSE)
    if(removeiRT) {
        input = .handleDecoyProteins(input, "StandardType", "iRT", FALSE)
        # TODO: maybe refactor later to match style of .handleSharedPeptides etc
    }
    input = .handleOxidationPeptides(input, "PeptideSequence", "+16", removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .handleFiltering(input, "Truncated", 0L, "smaller", "fill", NA_real_) # trick

    if (.checkDDA(input)) {
        input = .aggregateMonoisotopicPeaks(input)
        # TODO: maybe replace with .summarizeMultipleMeasurements (performance and logic). Can this be done later inside .summarize...?
    }
    
    input = .handleFiltering(input, "DetectionQValue", qvalue_cutoff, "smaller", "fill", 0.0, FALSE, filter_with_Qvalue)
    # TODO: conditionally check if DetectionQValue is in the input at the beginning
    feature_cols = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")
    input = .cleanByFeature(input, feature_cols, max, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature)
    input <- merge(input[, setdiff(colnames(input), c("Condition", "BioReplicate")), with = FALSE], annotation, by = "Run")
    input # TODO: change ProteinName to facor?
}
