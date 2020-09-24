#' Clean raw DIAUmpire files
#' @param msstats_object Object that inherits from MSstatsInputFiles class.
#' @param use_frag TRUE will use the selected fragment for each peptide. 'Selected_fragments' column is required.
#' @param use_pept TRUE will use the selected fragment for each protein 'Selected_peptides' column is required.
#' @return data.table
#' @keywords internal
.cleanRawDIAUmpire = function(msstats_object, use_frag, use_pept) {
    FragmentIon = ProteinName = PeptideSequence = . = ProteinKey = NULL
    frag_input = getInputFile(msstats_object, "Fragments")
    pept_input = getInputFile(msstats_object, "Peptides")
    prot_input = getInputFile(msstats_object, "Proteins")
    
    if (!is.element("Selected_fragments", colnames(pept_input)) | 
        !is.element("Selected_peptides", colnames(prot_input))) {
        msg = "Selected_fragments column is required. Please check it."
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    
    pept_cols = c("PeptideKey", "Proteins", "Selected_fragments")
    pept_input = pept_input[, pept_cols, with = FALSE]
    pept_input = pept_input[pept_input[["Proteins"]] != "", ]
    pept_input[["Selected_fragments"]] =  as.character(pept_input[["Selected_fragments"]])
    data.table::setnames(
        pept_input,
        colnames(pept_input), 
        c("PeptideSequence", "ProteinName", "FragmentIon"),
        skip_absent = TRUE)
    pept_input = pept_input[, lapply(.(FragmentIon), 
                                     function(x) unlist(data.table::tstrsplit(x, "|", fixed = TRUE))),
                            by = .(ProteinName, PeptideSequence)] 
    data.table::setnames(pept_input, "V1", "FragmentIon", skip_absent = TRUE)
    pept_input = pept_input[pept_input[["FragmentIon"]] != "", ]
    # TOOD: It's possible to extract a function here for prot/pept input  
    pept_input[["ProteinName"]] = as.character(pept_input[["ProteinName"]])
    
    prot_input = prot_input[, c("ProteinKey", "Selected_peptides"), with = FALSE]
    prot_input = prot_input[ProteinKey != "", ] 
    prot_input[["Selected_peptides"]] = as.character(prot_input[["Selected_peptides"]])
    prot_input[["Selected_peptides"]] = trimws(prot_input[["Selected_peptides"]],
                                               "both") # Is it needed?
    data.table::setnames(prot_input, colnames(prot_input),
                         c("ProteinName", "PeptideSequence"),
                         skip_absent = TRUE)
    prot_input = prot_input[, lapply(.(PeptideSequence),
                                     function(x) (unlist(data.table::tstrsplit(x, "|", fixed = TRUE)))),
                            by = .(ProteinName)]
    data.table::setnames(prot_input, "V1", "PeptideSequence", 
                         skip_absent = TRUE)
    prot_input = prot_input[PeptideSequence != "", ]
    
    prot_input[["ProteinName"]] = gsub(";", "", prot_input[["ProteinName"]])
    
    # TODO: move the comments to documentation
    if (use_frag & use_pept) {
        ## The number of peptides are the same in raw.pep3 or raw.pro3, but the number of proteins are different.
        ## usd these selected peptides, and use protein id in raw.pro3. because pro3 assign one protein id for each peptides
        ## now raw.pep3 can be used for final list but need to assign protein id from raw.pro3. It should be updated in raw.frag.
        ## raw.frag can have multiple protein ids.
        ## make the final list of selected proteins and peptides
        input = merge(pept_input[, c(FALSE, TRUE, TRUE), with = FALSE], 
                      prot_input, all.x = TRUE, by = "PeptideSequence",
                      sort = FALSE)
    } else if (use_frag & !use_pept) { 
        # only use selected fragment, and use all peptides, can use raw.pep3 only.
        ## but doublecheck protein id : protein id in raw.frag is always one id.
        ## need to find shared peptides
        input = pept_input
        input$ProteinName = gsub(";", "", input$ProteinName)
    } else {
        msg = "MSstats recommends to use at least selected fragments."
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    input[["FragmentIon"]] = gsub("\\+", "_", input[["FragmentIon"]])
    
    intensity_cols = grepl("Intensity", colnames(frag_input))
    frag_col_names = c("FragmentKey", "Protein", "Peptide", "Fragment")
    key_cols = colnames(frag_input) %in% frag_col_names
    frag_cols = intensity_cols | key_cols
    frag_input = frag_input[, frag_cols, with = FALSE]
    new_names = c("ProteinName", "PeptideSequence", "FragmentIon")
    data.table::setnames(frag_input, frag_col_names[2:4], new_names,
                         skip_absent = TRUE)
    frag_input = .fixColumnTypes(
        frag_input, character_columns = c(new_names, "FragmentKey"))
    
    if (!grepl("\\+", frag_input[["FragmentIon"]][1])) {
        # TODO: is always all or none true? Indicate in docs that + means charge
        last_char = nchar(frag_input[["FragmentKey"]])
        frag_input[["FragmentIon"]] = paste(
            frag_input[["FragmentIon"]],
            substr(frag_input[["FragmentKey"]], start = last_char, stop = last_char), 
            sep = "_")
    } else {
        frag_input[["FragmentIon"]] = gsub("\\+", "_", frag_input[["FragmentIon"]])
    }
    frag_input = frag_input[, (colnames(frag_input) != "FragmentKey"),
                            with = FALSE]
    
    setkey(input, ProteinName, PeptideSequence, FragmentIon)
    setkey(frag_input, ProteinName, PeptideSequence, FragmentIon)
    to_join = unique(frag_input[input, which = TRUE, allow.cartesian = TRUE, drop = FALSE])
    frag_input = frag_input[to_join]
    frag_input = melt(frag_input, 
                      id.vars = new_names, variable.name = "Run", 
                      value.name = "Intensity", value.factor = FALSE)
    frag_input[["Run"]] = gsub("_Intensity", "", frag_input[["Run"]])
    frag_input[!(is.na(ProteinName))]
}