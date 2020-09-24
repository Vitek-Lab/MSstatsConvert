#' Clean raw Progenesis output
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @param runs chr, vector of Run labels.
#' @param fix_colnames lgl, if TRUE, one of the rows will be used as colnames.
#' @return data.table
#' @keywords internal
.cleanRawProgenesis = function(msstats_object, runs, fix_colnames = TRUE) {
    Useinquantitation = NULL
    
    prog_input = getInputFile(msstats_object, "input")
    raw_abundance_col_id = which(colnames(prog_input) == "Rawabundance")
    if (is.element("Rawabundance", colnames(prog_input)) &
        is.element("Normalizedabundance", colnames(prog_input))) {
        norm_abundance_col_id = which(colnames(prog_input) == "Normalizedabundance")
        if ((raw_abundance_col_id - norm_abundance_col_id) != length(runs)) {
            msg = paste("** Please check annotation file. The numbers of MS runs", 
                        "in annotation and output are not matched.")
            getOption("MSstatsLog")("ERROR", msg)
            stop(msg)
        }
    }
    raw_abundances_col_ids = seq(raw_abundance_col_id,
                                 raw_abundance_col_id + length(runs) - 1)
    
    if (fix_colnames) {
        if (all(unique(prog_input[[1]][1:2]) == "")) {
            skip = 1:2
        } else {
            skip = 1
        }
        prog_input = prog_input[-skip, ]
        colnames(prog_input) = as.character(unlist(prog_input[1, ]))
        prog_input = prog_input[-1, ]
    }
    protein_col = .findAvailable(c("Protein", "Accession"), colnames(prog_input))
    cols = which(colnames(prog_input) %in% c(protein_col, "Modifications", 
                                             "Sequence", "Charge"))
    cols = c(cols, raw_abundances_col_ids)
    prog_input = prog_input[, cols, with = FALSE]
    colnames(prog_input) = .standardizeColnames(colnames(prog_input))
    data.table::setnames(prog_input, 
                         c(protein_col, "Charge"), 
                         c("ProteinName", "PrecursorCharge"),
                         skip_absent = TRUE)
    
    
    nonmissing_prot = !is.na(prog_input$ProteinName) & prog_input$ProteinName != ""
    nonmissing_pept = !is.na(prog_input$Sequence) & prog_input$Sequence != ""
    prog_input = prog_input[nonmissing_prot & nonmissing_pept, ]
    prog_input$PeptideSequence = paste(prog_input$Sequence,
                                       prog_input$Modifications,
                                       sep = "")
    if (is.element("Useinquantitation", colnames(prog_input))) {
        if (!is.logical(prog_input$Useinquantitation)) {
            prog_input$Useinquantitation = prog_input$Useinquantitation == "True"
        }
        prog_input = prog_input[(Useinquantitation), ]
        prog_input = prog_input[, !(colnames(prog_input) == "Useinquantitation"),
                                with = FALSE]
    }
    feature_cols = intersect(colnames(prog_input), 
                             c("ProteinName", "PeptideSequence",
                               "PrecursorCharge", "Fraction"))
    prog_cols = intersect(colnames(prog_input), 
                          c(feature_cols, as.character(runs)))
    prog_input = prog_input[, prog_cols, with = FALSE]
    prog_input = data.table::melt(prog_input,
                                  id.vars = feature_cols,
                                  measure.vars = runs,
                                  variable.name = "Run",
                                  value.name = "Intensity",
                                  value.factor = FALSE)
    prog_input$Intensity = as.numeric(as.character(prog_input$Intensity))
    prog_input
}
