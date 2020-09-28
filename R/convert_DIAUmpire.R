#' Clean raw DIAUmpire files
#' @param msstats_object Object that inherits from MSstatsInputFiles class.
#' @param use_frag TRUE will use the selected fragment for each peptide. 
#' 'Selected_fragments' column is required.
#' @param use_pept TRUE will use the selected fragment for each protein 
#' 'Selected_peptides' column is required.
#' @return data.table
#' @keywords internal
.cleanRawDIAUmpire = function(msstats_object, use_frag, use_pept) {
    FragmentIon = ProteinName = PeptideSequence = Selected_peptides = NULL
    Selected_fragments = Proteins = ProteinKey = FragmentKey = Run = NULL
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
    pept_input = pept_input[Proteins != "", ]
    pept_input[, Selected_fragments := as.character(Selected_fragments)]
    data.table::setnames(
        pept_input,
        colnames(pept_input), 
        c("PeptideSequence", "ProteinName", "FragmentIon"),
        skip_absent = TRUE)
    pept_input = pept_input[
        , lapply(list(FragmentIon), 
                 function(x) unlist(data.table::tstrsplit(x, "|", 
                                                          fixed = TRUE))),
        by = c("ProteinName", "PeptideSequence")
        ] 
    data.table::setnames(pept_input, "V1", "FragmentIon", skip_absent = TRUE)
    pept_input = pept_input[pept_input[["FragmentIon"]] != "", ]
    pept_input[["ProteinName"]] = as.character(pept_input[["ProteinName"]])
    
    prot_input = prot_input[, c("ProteinKey", "Selected_peptides"), 
                            with = FALSE]
    prot_input = prot_input[ProteinKey != "", ] 
    prot_input[, Selected_peptides := as.character(Selected_peptides)]
    prot_input[, Selected_peptides := trimws(Selected_peptides, "both")]
    data.table::setnames(prot_input, colnames(prot_input),
                         c("ProteinName", "PeptideSequence"),
                         skip_absent = TRUE)
    prot_input = prot_input[
        , lapply(list(PeptideSequence),
                 function(x) (unlist(data.table::tstrsplit(x, "|", 
                                                           fixed = TRUE)))),
        by = c("ProteinName")]
    data.table::setnames(prot_input, "V1", "PeptideSequence", 
                         skip_absent = TRUE)
    prot_input = prot_input[PeptideSequence != "", ]
    prot_input[, ProteinName := gsub(";", "", ProteinName)]

    if (use_frag & use_pept) {
        input = merge(pept_input[, c(FALSE, TRUE, TRUE), with = FALSE], 
                      prot_input, all.x = TRUE, by = "PeptideSequence",
                      sort = FALSE)
    } else if (use_frag & !use_pept) { 
        input = pept_input
        input$ProteinName = gsub(";", "", input$ProteinName)
    } else {
        msg = "MSstats recommends to use at least selected fragments."
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    input[, FragmentIon := gsub("\\+", "_", FragmentIon)]

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
        last_char = nchar(frag_input[["FragmentKey"]])
        frag_input[, FragmentIon := paste(FragmentIon, 
                                          substr(FragmentKey, last_char, 
                                                 last_char), sep = "_")]
    } else {
        frag_input[, FragmentIon := gsub("\\+", "_", FragmentIon)]
    }
    frag_input = frag_input[, (colnames(frag_input) != "FragmentKey"),
                            with = FALSE]
    
    setkey(input, ProteinName, PeptideSequence, FragmentIon)
    setkey(frag_input, ProteinName, PeptideSequence, FragmentIon)
    to_join = unique(frag_input[input, which = TRUE, 
                                allow.cartesian = TRUE, drop = FALSE])
    frag_input = frag_input[to_join]
    frag_input = melt(frag_input, 
                      id.vars = new_names, variable.name = "Run", 
                      value.name = "Intensity", value.factor = FALSE)
    frag_input[, Run := gsub("_Intensity", "", Run)]
    frag_input[!(is.na(ProteinName))]
}