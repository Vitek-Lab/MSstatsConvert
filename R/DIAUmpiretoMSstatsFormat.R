#' Import DIA-Umpire files 
#' 
#' @inheritParams .documentFunction
#' @param raw.frag name of FragSummary_date.xls data, which includes feature-level data.
#' @param raw.pep name of PeptideSummary_date.xls data, which includes selected fragments information.
#' @param raw.pro name of ProteinSummary_date.xls data, which includes selected peptides information.
#' @param annotation name of annotation data which includes Condition, BioReplicate, Run information.
#' @param useSelectedFrag TRUE will use the selected fragment for each peptide. 'Selected_fragments' column is required.
#' @param useSelectedPep TRUE will use the selected peptide for each protein. 'Selected_peptides' column is required.
#' 
#' @return data.frame with the required format of MSstats.
#'
#' @author Meena Choi, Olga Vitek 
#'
#' @export
#' 

DIAUmpiretoMSstatsFormat <- function(
    raw.frag, raw.pep, raw.pro, annotation, useSelectedFrag = TRUE,
    useSelectedPep = TRUE, fewMeasurements = "remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max){
    
    fewMeasurements = .isLegalValue(fewMeasurements, c("keep", "remove"))
    # TODO: annotation checks
    
    input = .cleanRawDIAUmpire(raw.frag, raw.pep, raw.pro, useSelectedFrag,
                               useSelectedPep)
    input = .handleSharedPeptides(input, TRUE) # this function always removes and does it earlier than others
    # TODO: can I move removing shared peptides here?
    feature_cols = c("PeptideSequence", "FragmentIon")
    input = .cleanByFeature(input, feature_cols, summaryforMultipleRows, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Feature)
    input = merge(input, annotation, by = "Run")
    
    input = .fillValues(input, c("PrecursorCharge" = NA
                                 "ProductCharge" = NA
                                 "IsotopeLabelType" = "L"))
    input # ProteinName as factor?
}
