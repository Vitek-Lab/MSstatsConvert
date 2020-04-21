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

DIAUmpiretoMSstatsFormat = function(
    raw.frag, raw.pep, raw.pro, annotation, useSelectedFrag = TRUE,
    useSelectedPep = TRUE, fewMeasurements = "remove",
    removeProtein_with1Feature = FALSE, summaryforMultipleRows = max,
    use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
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
    
    input = .fillValues(input, c("PrecursorCharge" = NA,
                                 "ProductCharge" = NA,
                                 "IsotopeLabelType" = "L"))
    input # ProteinName as factor?
}


.cleanRawDIAUmpire = function(frag_input, pept_input, prot_input,
                              use_frag, use_pept) {
    if(!is.element("Selected_fragments", pept_input)) {
        stop("** Selected_fragments column is required. Please check it.")
    }
    if(!is.element("Selected_peptides", prot_input)) {
        stop("** Selected_peptides column is required. Please check it.")
    }
    
    pept_input = data.table::as.data.table(pept_input)
    pept_cols = c("Peptide.Key", "Proteins", "Selected_fragments")
    pept_input = pept_input[, pept_cols, with = FALSE]
    pept_input = pept_input[pept_input[["Proteins"]] != "", ]
    pept_input[["Selected_fragments"]] =  as.character(pept_input[["Selected_fragments"]])
    colnames(pept_input) = .updateColnames(
        pept_input,
        colnames(pept_input), 
        c("PeptideSequence", "ProteinName", "FragmentIon"))
    pept_input = pept_input[, lapply(.(FragmentIon), 
                                     function(x) unlist(tstrsplit(x, "|", fixed = TRUE))),
                            by = .(ProteinName, PeptideSequence)] 
    # TODO: test this data.table idiom properly!
    colnames(pept_input) = .updateColnames(pept_input, "V1", "FragmentIon")
    pept_input = pept_input[pept_input[["FragmentIon"]] != "", ]
    # TOOD: It's possible to extract a function here for prot/pept input  
    pept_input[["ProteinName"]] = as.character(pept_input[["ProteinName"]])
    
    prot_input = data.table::as.data.table(prot_input)
    prot_input = prot_input[, c("Protein.Key", "Selected_peptides"), with = FALSE]
    prot_input = prot_input[prot_input[["Protein.Key"]] != "", ] # change to the decoy function?
    prot_input[["Selected_peptides"]] = as.character(prot_input[["Selected_peptides"]])
    # prot_input[["Selected_peptides"]] = trimws(prot_input[["Selected_peptides"]],
    #                                            "both") # Is it needed?
    colnames(prot_input) = .updateColnames(prot_input, colnames(prot_input),
                                           c("ProteinName", "PeptideSequence"))
    prot_input = prot_input[, lapply(.(PeptideSequence),
                                     function(x) (unlist(tstrsplit(x, "|", fixed = TRUE)))),
                            by = .(ProteinName)]
    colnames(prot_input) = .updateColnames(prot_input, "V1", "PeptideSequence")
    prot_input = prot_input[prot_input[["PeptideSequence"]] != "", ]
    prot_input[["ProteinName"]] = gsub(";", "", prot_input[["ProteinName"]])
    
    # TODO: move the comments to documentation
    if (use_frag & use_pept) {
        ## The number of peptides are the same in raw.pep3 or raw.pro3, but the number of proteins are different.
        ## usd these selected peptides, and use protein id in raw.pro3. because pro3 assign one protein id for each peptides
        ## now raw.pep3 can be used for final list but need to assign protein id from raw.pro3. It should be updated in raw.frag.
        ## raw.frag can have multiple protein ids.
        ## make the final list of selected proteins and peptides
        input = merge(pept_input[, c(FALSE, TRUE, TRUE), with = FALSE], 
                      prot_input, all.x = TRUE, by = "PeptideSequence")
        input[["ProteinName"]] = as.character(input[["ProteinName"]])
    } else if (use_frag & !use_pept) { 
        # only use selected fragment, and use all peptides, can use raw.pep3 only.
        ## but doublecheck protein id : protein id in raw.frag is always one id.
        ## need to find shared peptides
        input = pept_input
    } else if (!use_frag & !use_pept) {
        stop('** MSstats recommends to use at least selected fragments.')
    }
    input[["FragmentIon"]] = gsub("\\+", "_", input[["FragmentIon"]])
    
    frag_input = data.table::as.data.table(frag_input)
    frag_input[["Protein"]] = as.character(frag_input[["Protein"]])
    intensity_cols = grepl("Intensity", colnames(frag_input))
    frag_col_names = c("Fragment.Key", "Protein", "Peptide", "Fragment")
    key_cols = colnames(frag_input) %in% frag_col_names
    frag_cols = intensity_cols | key_cols
    frag_input = frag_input[, frag_cols, with = FALSE]
    new_names = c("ProteinName", "PeptideSequence", "FragmentIon")
    colnames(frag_input) = .updateColnames(
        frag_input, frag_col_names[2:4], new_names)
    frag_input = .fixColumnTypes(
        frag_input, character_columns = c(new_names, "Fragment.Key"))
    
    if (!grepl("\\+", frag_input[["FragmentIon"]][1])) {
        # TODO: is alway all or none true? Indicate in docs that + means charge
        last_char = nchar(frag_input[["Fragment.Key"]])
        frag_input[["FragmentIon"]] = paste(
            frag_input[["FragmentIon"]],
            substr(frag_input[["Fragment.Key"]], start = last_char, stop = last_char), 
            sep = "_")
    } else {
        frag_input[["FragmentIon"]] = gsub("\\+", "_", frag_input[["FragmentIon"]])
    }
    frag_input = frag_input[, (colnames(frag_input) != "Fragment.Key"),
                            with = FALSE]
    
    setkey(input, ProteinName, PeptideSequence, FragmentIon)
    setkey(frag_input, ProteinName, PeptideSequence, FragmentIon)
    to_join = unique(frag_input[input, which = TRUE, allow.cartesian = TRUE])
    frag_input = frag_input[to_join]
    frag_input = melt(frag_input, 
                      id.vars = new_names,
                      variable.name = "Run", value.name = "Intensity")
    frag_input[["Run"]] = gsub("_Intensity", "", frag_input[["Run"]])
    frag_input 
}