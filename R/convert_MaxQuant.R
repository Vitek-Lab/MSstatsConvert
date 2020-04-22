#' Import MaxQuant files
#' 
#' @inheritParams .documentFunction
#' @inheritParams .setMS
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

MaxQtoMSstatsFormat = function(
    evidence, annotation, proteinGroups, proteinID = "Proteins", 
    useUniquePeptide = TRUE, summaryforMultipleRows = max, 
    fewMeasurements = "remove", removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
    
    .setMSstatsLogger(use_log_file, append, verbose)
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
    
    input = .getDataTable(evidence)
    proteinGroups = .getDataTable(proteinGroups)
    input = .cleanRawMaxQuant(input, proteinGroups, proteinID)
    input = .handleOxidationPeptides(input, "PeptideSequence", 
                                     "M", removeMpeptides)
    input = .handleOxidationPeptides(input, "Modifications", "Oxidation",
                                     removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .cleanByFeature(input, c("PeptideSequence", "PrecursorCharge"), summaryforMultipleRows, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide,
                                           c("PeptideSequence", "PrecursorCharge"))
    input = merge(input, annotation, by = "Run") # by , "IsotopeLabelType?
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 "IsotopeLabelType"  =  "L"))
    new("MSstatsValidated", as.data.frame(input)) # Convert ProteinName and PeptideSequence to factor?
}


.cleanRawMaxQuant = function(mq_input, mq_pg, protein_id_col) {
    colnames(mq_input) = .standardizeColnames(colnames(mq_input))
    colnames(mq_pg) = .standardizeColnames(colnames(mq_pg))
    filter_cols = c("Contaminant", "Potential.contaminant", 
                    "Reverse", "Only.identified.by.site")
    mq_input = .filterManyColumns(mq_input, filter_cols, "+")
    mq_pg = .filterManyColumns(mq_pg, filter_cols, "+")
    mq_input[["Protein.group.IDs"]] = suppressWarnings(as.integer(as.character(mq_input[["Protein.group.IDs"]])))
    # Need to check if 'id' in proteinGroups.txt and 'Protein.group.IDs' in evidence are the same. 
    # Possible to have some combination in Protein.group.IDs in infile, such as 64;1274;1155;1273 instead of 64, 1274.. separately. 
    # Combination of some ids seems not to be used for intensity
    # TODO: move this comment to docs.
    ## then take proteins which are included
    mq_input = mq_input[mq_input[["Protein.group.IDs"]] %in% unique(mq_pg[["id"]]), ]
    ## then use 'protein.IDs' in proteinGroups.txt
    ## because if two 'proteins' in evidence.txt are used in one protein ID, need to use certain protein name in infile.
    ## for example, protein.IDs in proteinGroups.txt are P05204;O00479. but, two 'proteins in evidence.txt, such as P05204;O00479, and P05204.
    # TODO: move this comment to docs
    protein_names = unique(mq_pg[,c("Protein.IDs", "id")])
    colnames(protein_names) = c("uniqueProteins", "Protein.group.IDs")
    mq_input = merge(mq_input, protein_names, by = "Protein.group.IDs")
    
    mq_cols = c("Sequence", "Modified.sequence", 
                "Modifications", "Charge", "Run", "Intensity", 
                "Fraction", "TechReplicate", "Raw.file") 
    # TODO: + Retention.time? - Modifications?
    protein_id = ifelse(protein_id_col == "Proteins", "uniqueProteins",
                        "Leading.razor.protein")
    mq_cols = intersect(c(mq_cols, protein_id),
                        colnames(mq_input))
    mq_input = mq_input[, mq_cols, with = FALSE]
    mq_input[["Modified.sequence"]] = gsub("_", "", mq_input[["Modified.sequence"]])
    mq_input = mq_input[!is.na(mq_input[["Intensity"]]), ]
    colnames(mq_input) = .updateColnames(
        mq_input, 
        c(protein_id, "Modified.sequence", "Charge", "Raw.file"), 
        c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run"))
    mq_input
}
