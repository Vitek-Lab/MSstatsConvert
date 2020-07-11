#' Clean raw Progenesis output
#' @param msstats_object an object of class `MSstatsSpectroMineFiles`.
#' @param runs chr, vector of Run labels.
#' @param fix_colnames lgl, if TRUE, one of the rows will be used as colnames.
#' @return data.table
#' @keywords internal
.cleanRawProgenesis = function(msstats_object, runs, fix_colnames = TRUE) {
    Useinquantitation = NULL
    
    prog_input = getInputFile(msstats_object, "input")
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
    colnames(prog_input) = .standardizeColnames(colnames(prog_input))
    protein_col = .findAvailable(c("Protein", "Accession"), 
                                 colnames(prog_input))
    colnames(prog_input) = .updateColnames(prog_input, 
                                           c(protein_col, "Charge"), 
                                           c("ProteinName", "PrecursorCharge"))
    
    nonmissing_prot = !is.na(prog_input$ProteinName) & prog_input$ProteinName != ""
    nonmissing_pept = !is.na(prog_input$Sequence) & prog_input$Sequence != ""
    prog_input = prog_input[nonmissing_prot & nonmissing_pept, ]
    prog_input$PeptideSequence = paste(prog_input$Sequence,
                                       prog_input$Modifications,
                                       sep = "")
    prog_input = prog_input[!duplicated(prog_input), ] # dubious performance-wise
    if (is.element("Useinquantitation", colnames(prog_input))) {
        if (!is.logical(prog_input$Useinquantitation)) {
            prog_input$Useinquantitation = prog_input$Useinquantitation == "True"
        }
        prog_input = prog_input[(Useinquantitation), ]
        # TODO: consider character version if these files ever import this columns as character
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
    prog_input$Intensity = as.numeric(prog_input$Intensity)
    prog_input
}
