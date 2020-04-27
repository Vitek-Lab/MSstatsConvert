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

MaxQtoMSstatsFormat = function(
    evidence, annotation, proteinGroups, proteinID = "Proteins", 
    useUniquePeptide = TRUE, summaryforMultipleRows = max, 
    fewMeasurements = "remove", removeMpeptides = FALSE,
    removeOxidationMpeptides = FALSE, removeProtein_with1Peptide = FALSE,
    use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
    .setMSstatsLogger(use_log_file, append, verbose)
    # .checkConverterParams()
    
    input = .cleanRawMaxQuant(evidence, proteinGroups, proteinID)
    input = .makeAnnotation(input, .getDataTable(annotation), "Run" = "Rawfile")
        
    input = .handleOxidationPeptides(input, "PeptideSequence", 
                                     "M", removeMpeptides)
    input = .handleOxidationPeptides(input, "Modifications", "Oxidation",
                                     removeOxidationMpeptides)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .cleanByFeature(input, c("PeptideSequence", "PrecursorCharge"), summaryforMultipleRows, fewMeasurements)
    input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide,
                                           c("PeptideSequence", "PrecursorCharge"))
    input = .mergeAnnotation(input, annotation) 
    input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                                 "IsotopeLabelType"  =  "L"))
    input
}


#' Clean raw output from MaxQuant
#' @param mq_input evidence file from MaxQuant or a path to it.
#' @param mq_pg proteinGroups file from MaxQuant or a path to it.
#' @param protein_id_col chr, name of a column with names of proteins.
#' @param ... optional, additional parameters for data.table::fread
#' @return data.table
#' @keywords internal
.cleanRawMaxQuant = function(mq_input, mq_pg, protein_id_col, ...) {
    mq_input = .getDataTable(mq_input, ...)
    mq_pg = .getDataTable(mq_pg, ...)
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


MaxQtoMSstatsTMTFormat = function(
    evidence, proteinGroups, annotation, which.proteinid = 'Proteins',
    rmProt_Only.identified.by.site = FALSE, useUniquePeptide = TRUE,
    rmPSM_withMissing_withinRun = FALSE, rmPSM_withfewMea_withinRun = TRUE,
    rmProtein_with1Feature = FALSE, summaryforMultipleRows = sum,
    use_log_file = TRUE, append = TRUE, verbose = TRUE
) {
    .setMSstatsLogger(use_log_file, append, verbose)
    # .checkConverterParams()

    input = .cleanRawMaxQuantTMT(evidence)
    annotation = .makeAnnotation(input, .getDataTable(annotation))

    feature_cols = c("PeptideSequence", "Charge")
    input = .removeMissingAllChannels(input)
    input = .handleSharedPeptides(input, useUniquePeptide)
    input = .cleanByFeatureTMT(input, feature_cols, summaryforMultipleRows, 
                               rmPSM_withfewMea_withinRun, rmPSM_withMissing_withinRun)
    input = .mergeAnnotation(input, annotation)
    input = .handleSingleFeaturePerProtein(input, rmProtein_with1Feature, "PSM")
    input = .handleFractions(input, annotation)
    input[, c("ProteinName", "PeptideSequence", "Charge", "PSM", "Mixture", 
              "TechRepMixture", "Run", "Channel", "BioReplicate", "Condition", "Intensity")]
    
}

.cleanRawMaxQuantTMT = function(mq_input, mq_pg, remove_by_site = FALSE,
                                channel_columns = "Reporter.intensity.corrected",
                                ...) {
    mq_input = .getDataTable(mq_input, ...)
    colnames(mq_input) = .standardizeColnames(colnames(mq_input))
    mq_pg = .getDataTable(mq_pg, ...)
    colnames(mq_pg) = .standardizeColnames(colnames(mq_pg))
    
    mq_input = .filterManyColumns(mq_input, 
                                  c("Contaminant", "Potential.contaminant", "Reverse"), 
                                  "+")
    msg = "+ Contaminant, + Reverse, + Only.identified.by.site, proteins are removed."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    
    mq_input = .filterExact(mq_input, "Only.identified.by.site", "+",
                            is.element("Only.identified.by.site", colnames(mq_input)), 
                            remove_by_site)
    msg = "+ Only.identified.by.site are removed."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg) # TODO: make this conditional
    
    mq_pg = .filterManyColumns(mq_pg, 
                               c("Contaminant", "Potential.contaminant", "Reverse"), 
                               "+")
    mq_pg = .filterExact(mq_pg, "Only.identified.by.site", 
                         is.element("Only.identified.by.site", 
                                    colnames(mq_pg)), 
                         remove_by_site)
    mq_input$Protein.group.IDs = as.numeric(mq_input$Protein.group.IDs)
    mq_input = merge(mq_input, 
                     unique(mq_pg[, .(uniquefromProteinGroups = Protein.IDs,
                                      Protein.group.IDs = id)]),
                     by = "Protein.group.IDs") # is unique() necessary? + check this code
    protein_ID = .findAvailable(c("Proteins", "Leading.proteins", 
                                  "Leading.razor.protein", "Gene.names"),
                                colnames(mq_input), "Proteins")
    channels = .getChannelColumns(colnames(mq_input), channel_columns)
    mq_input = mq_input[, c(protein_ID, "Modified.sequence", "Charge", "Raw.file",
                            "Score", channels), with = FALSE]
    colnames(mq_input) = .updateColnames(mq_input, 
                                         c(protein_ID, "Modified.sequence", "Raw.file"),
                                         c("ProteinName", "PeptideSequence", "Run"))
    mq_input[["PeptideSequence"]] = gsub("_", "", mq_input[["PeptideSequence"]])
    mq_input$PSM = paste(mq_input$PeptideSequence, 
                         mq_input$PeptideCharge, 
                         1:nrow(mq_input),
                         sep = "_")
    mq_input = melt(mq_input, measure.vars = channels,
                    id.vars = c("ProteinName", "PeptideSequence", "Charge", "PSM", "Run", "Score"),
                    variable.name = "Channel", value.name = "Intensity")
    mq_input$Channel = gsub(channel_columns, "channel", mq_input$Channel)
    mq_input
}
