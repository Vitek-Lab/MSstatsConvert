#' Import MaxQuant files
#' 
#' @inheritParams .documentFunction
#' @param evidence name of 'evidence.txt' data, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Raw.file, Condition, BioReplicate, Run, IsotopeLabelType information.
#' @param proteinGroups name of 'proteinGroups.txt' data. It needs to matching protein group ID. If proteinGroups=NULL, use 'Proteins' column in 'evidence.txt'.
#' @param proteinID 'Proteins'(default) or 'Leading.razor.protein' for Protein ID.
#' 
#' @return data.frame with the required format of MSstats.
#' 
#' @note Warning: MSstats does not support for metabolic labeling or iTRAQ experiments.
#' 
#' @author Meena Choi, Olga Vitek. 
#' 
#' @export
#' 

MaxQtoMSstatsFormat <- function(
    evidence, annotation, proteinGroups, proteinID = "Proteins", 
    useUniquePeptide = TRUE, summaryforMultipleRows = max, 
    fewMeasurements = "remove", removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE) {
    
    fewMeasurements = .isLegalValue(fewMeasurements, 
                                    legal_values = c("remove", "keep"))
    proteinID = .isLegalValue(proteinID, 
                              legal_values = c("Proteins", 
                                               "Leading.razor.protein"))
    annotation = .makeAnnotation(
        annotation,
        c("Raw.file" = "Run", "Condition" = "Condition", "BioReplicate" = "BioReplicate",
          "IsotopeLabelType" = "IsotopeLabelType")
    )

    input = .cleanRawMaxQuant(evidence, proteinGroups)
    input = .handleOxidationPeptides(input, "Modified.sequence", 
                                     "M", removeMpeptides)
    input = .handleOxidationPeptides(input, "Modifications", "Oxidation",
                                     removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .cleanByFeature(input, c("PeptideSequence", "PrecursorCharge"), summaryforMultipleRows, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide)
    input = merge(input, annotation, by = c("Run", "IsotopeLabelType"))
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 IsotopeLabelType  =  "L"))
    input # Convert ProteinName and PeptideSequence to factor?
    # TODO: check with actual annotation from MaxQuant
}
