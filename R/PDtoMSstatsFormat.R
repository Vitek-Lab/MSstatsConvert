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
    
    fewMeasurements = .isLegalValue(fewMeasurements, 
                                    legal_values = c("remove", "keep"))
    
    annotation = .makeAnnotation(annotation, 
                                 c("Run" = "Run", "Condition" = "Condition", 
                                   "BioReplicate" = "BioReplicate"))
    input = .cleanRawPD(input, which.proteinid, which.quantification, 
                        which.sequence, useNumProteinsColumn)
    feature_cols = c("PeptideModifiedSequence", "PrecursorCharge", "FragmentIon")
    input = .handleOxidationPeptides(input, "PeptideModifiedSequence", 
                                      "Oxidation", removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide,
                                  peptide_column = "PeptideModifiedSequence")
    input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                            fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, feature_cols, 
                                           removeProtein_with1Peptide)
    input = merge(input, annotation, by = "Run")
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 "IsotopeLabelType" = "L"))
    input
}
