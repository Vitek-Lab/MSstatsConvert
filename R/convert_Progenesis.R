#' Import Progenesis files
#' 
#' @inheritParams .documentFunction
#' @param input name of Progenesis output, which is wide-format. 'Accession', 'Sequence', 'Modification', 'Charge' and one column for each run are required.
#' @param annotation name of 'annotation.txt' or 'annotation.csv' data which includes Condition, BioReplicate, Run information. It will be matched with the column name of input for MS runs.
#'
#' @return data.frame with the required format of MSstats.
#' 
#' @author Meena Choi, Olga Vitek, Ulrich Omasits
#' 
#' @export
#' 

ProgenesistoMSstatsFormat = function(
  input, annotation, useUniquePeptide = TRUE, summaryforMultipleRows = max,
  fewMeasurements = "remove", removeOxidationMpeptides = FALSE, 
  removeProtein_with1Peptide = FALSE,
  use_log_file = TRUE, append = FALSE, verbose = TRUE
) {
  fewMeasurements = .isLegalValue(fewMeasurements, c("remove", "keep"))
  # Check go here
  annotation = .makeAnnotation(annotation, 
                               c("Run" = "Run", "Condition" = "Condition",
                                 "BioReplicate" = "BioReplicate"))
  
  input = .cleanRawProgenesis(input, runs)
  input = .handleOxidationPeptides(input, "PeptideModifiedSequence", 
                                   "Oxidation", removeOxidationMpeptides)
  input = .handleSharedPeptides(input, useUniquePeptide,
                                peptide_column = "PeptideModifiedSequence")
  feature_cols = c("PeptideModifiedSequence", "PrecursorCharge")
  input = .cleanByFeature(input, feature_cols, summaryforMultipleRows,
                          fewMeasurements)
  input = .handleSingleFeaturePerProtein(input, removeProtein_with1Peptide,
                                         feature_cols)
  input = merge(input, annotation, by = "Run")
  input = .fillValues(input, c("FragmentIon" = NA, "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
  input
}


#' Clean raw Progenesis output
#' @param prog_input Progenesis report or a path to it.
#' @param runs chr, vector of Run labels.
#' @param fix_colnames lgl, if TRUE, one of the rows will be used as colnames.
#' @param ... optional, additional parameters to data.table::fread.
#' @return data.table
#' @keywords internal
.cleanRawProgenesis = function(prog_input, runs, fix_colnames = TRUE, ...) {
  prog_input = .getDataTable(prog_input, ...)
  if (fix_colnames) {
    prog_input = prog_input[-(1:2), ]
    colnames(prog_input) = unlist(prog_input[1, ])
    prog_input = prog_input[-1, ]
  }
  colnames(prog_input) = .standardizeColnames(colnames(prog_input))
  protein_col = .findAvailable(c("Protein", "Accession"), 
                               colnames(prog_input))
  colnames(prog_input) = .updateColnames(prog_input, 
                                         c(protein_col, "Charge"), 
                                         c("ProteinName", "PrecursorCharge"))
  
  nonmissing_prot = !is.na(prog_input[["ProteinName"]]) & prog_input[["ProteinName"]] != ""
  nonmissing_pept = !is.na(prog_input[["Sequence"]]) & prog_input[["Sequence"]] != ""
  prog_input = prog_input[nonmissing_prot & nonmissing_pept, ]
  prog_input[["PeptideModifiedSequence"]] = paste(prog_input[["Sequence"]],
                                                  prog_input[["Modifications"]],
                                                  sep = "")
  prog_input = prog_input[!duplicated(prog_input), ] # dubious performance-wise
  if (is.element("Use.in.quantitation", colnames(prog_input))) {
    if (!is.logical(prog_input[["Use.in.quantitation"]])) {
      prog_input[["Use.in.quantitation"]] = prog_input[["Use.in.quantitation"]] == "True"
    }
    prog_input = prog_input[prog_input$Use.in.quantitation, ]
    # TODO: consider character version if these files ever import this columns as character
    prog_input = prog_input[, !(colnames(prog_input) == "Use.in.quantitation"),
                            with = FALSE]
  }
  feature_cols = intersect(colnames(prog_input), 
                           c("ProteinName", "PeptideModifiedSequence",
                             "PrecursorCharge", "Fraction"))
  prog_cols = intersect(colnames(prog_input), c(feature_cols, runs))
  prog_input = prog_input[, prog_cols, with = FALSE]
  prog_input = data.table::melt(prog_input,
                                id.vars = feature_cols,
                                measure.vars = runs,
                                variable.name = "Run",
                                value.name = "Intensity",
                                value.factor = FALSE)
  prog_input[["Intensity"]] = as.numeric(prog_input[["Intensity"]])
  prog_input
}
