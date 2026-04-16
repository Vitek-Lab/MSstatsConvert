#' Clean raw Spectronaut output.
#' @param msstats_object an object of class `MSstatsSpectronautFiles`.
#' @inheritParams SpectronauttoMSstatsFormat
#' @return `data.table`
#' @keywords internal
.cleanRawSpectronaut = function(msstats_object, intensity,
                                calculateAnomalyScores,
                                anomalyModelFeatures,
                                peptideSequenceColumn = "EG.ModifiedSequence",
                                heavyLabels = NULL) {
  FFrgLossType = FExcludedFromQuantification = NULL

  spec_input = getInputFile(msstats_object, "input")
  .validateSpectronautInput(spec_input, peptideSequenceColumn)
  spec_input = .addSpectronautColumnsIfMissing(spec_input)
  spec_input = spec_input[FFrgLossType == "noloss", ]
  
  f_charge_col = .findAvailable(c("FCharge", "FFrgZ"), colnames(spec_input), fall_back = "FCharge")
  pg_qval_col = .findAvailable(c("PGQvalue"), colnames(spec_input))
  interference_col = .findAvailable(c("FPossibleInterference"),
                                    colnames(spec_input))
  exclude_col = .findAvailable(c("FExcludedFromQuantification"),
                               colnames(spec_input))
  intensity_col = .resolveSpectronautIntensityColumn(intensity, colnames(spec_input))
  peptide_col = .findAvailable(.standardizeColnames(peptideSequenceColumn),
                               colnames(spec_input))
  protein_col = .findAvailable(c("PGProteinGroups", "PGProteinAccessions"),
                               colnames(spec_input), fall_back = "PGProteinGroups")

  cols = c(protein_col, peptide_col, "FGCharge", "FFrgIon",
           f_charge_col, "RFileName", "RCondition", "RReplicate",
           "EGQvalue", pg_qval_col, interference_col, exclude_col,
           intensity_col)
  if (calculateAnomalyScores){
    cols = c(cols, anomalyModelFeatures)
  }
  cols = intersect(cols, colnames(spec_input))
  spec_input = spec_input[, cols, with = FALSE]

  data.table::setnames(
    spec_input,
    c(protein_col, peptide_col, "FGCharge", "FFrgIon",
      f_charge_col, "RFileName", intensity_col,
      "RCondition", "RReplicate"),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Intensity", "Condition", "BioReplicate"),
    skip_absent = TRUE)

  spec_input = .assignSpectronautIsotopeLabelType(
    spec_input, heavyLabels)

  .logSuccess("Spectronaut", "clean")
  spec_input
}

#' Helper method to validate input has necessary columns
#' @param spec_input dataframe input
#' @param peptideSequenceColumn character, name of the column containing peptide 
#' sequences, passed from user
#' @noRd
.validateSpectronautInput = function(spec_input, peptideSequenceColumn) {
    # Only FGCharge is truly required
    required_columns = c("FGCharge")
    missing_columns = setdiff(required_columns, colnames(spec_input))
    if (length(missing_columns) > 0) {
        msg = paste("The following columns are missing from the input data:",
                    paste(missing_columns, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    # Ensure at least one protein name column is present
    if (!any(c("PGProteinGroups", "PGProteinAccessions") %in% colnames(spec_input))) {
        msg = paste("The following columns are missing from the input data:",
                    "PGProteinGroups")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    # Ensure at least one protein name column is present
    if (!(.standardizeColnames(peptideSequenceColumn) %in% colnames(spec_input))) {
        msg = paste("The following column are missing from the input data:",
                    peptideSequenceColumn)
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
}

#' Add synthetic columns that are absent in protein turnover Spectronaut reports.
#'
#' Spectronaut's protein turnover export omits several columns that the
#' standard MSstats pipeline expects.  Rather than requiring callers to patch
#' the data frame before passing it in (the "hacks" documented in the design
#' document), we synthesize sensible defaults here so the rest of the pipeline
#' sees a consistent schema.
#'
#' Columns synthesized when absent:
#' \describe{
#'   \item{FFrgLossType}{"noloss" — the filter \code{spec_input[FFrgLossType == "noloss"]}
#'     keeps all rows, which is the correct behavior when the column is missing.}
#'   \item{FExcludedFromQuantification}{FALSE — no rows are excluded.}
#'   \item{FFrgIon}{NA — fragment ion identity is not available at the precursor
#'     level used by MS1-based protein turnover workflows.}
#'   \item{FCharge}{NA — product charge is not applicable when fragment-level
#'     columns are absent.}
#' }
#'
#' @param spec_input `data.table` with standardized column names.
#' @return `data.table` with missing columns added.
#' @keywords internal
#' @noRd
.addSpectronautColumnsIfMissing = function(spec_input) {
  if (!("FFrgLossType" %in% colnames(spec_input))) {
    spec_input[, FFrgLossType := "noloss"]
  }
  if (!("FExcludedFromQuantification" %in% colnames(spec_input))) {
    spec_input[, FExcludedFromQuantification := FALSE]
  }
  if (!("FFrgIon" %in% colnames(spec_input))) {
    spec_input[, FFrgIon := NA_character_]
  }
  if (!("FCharge" %in% colnames(spec_input))) {
    spec_input[, FCharge := NA_integer_]
  }
  spec_input
}


#' Resolve the Spectronaut intensity column from user input.
#'
#' Accepts either a legacy enum alias or a raw (standardized) column name.
#' The legacy aliases map to their canonical standardized column names so that
#' old code continues to work unchanged.  When the user passes a raw column
#' name (e.g. \code{"FGMS1Quantity"} or \code{"FGMS2Quantity"}), it is used
#' directly after verifying that the column exists.
#'
#' @param intensity Character scalar: enum alias or raw column name.
#' @param available_cols Character vector of available standardized column names.
#' @return The resolved standardized column name.
#' @keywords internal
#' @noRd
.resolveSpectronautIntensityColumn = function(intensity, available_cols) {
  legacy_mapping = c(
    "PeakArea"           = "FPeakArea",
    "NormalizedPeakArea" = "FNormalizedPeakArea"
  )

  if (intensity %in% names(legacy_mapping)) {
    resolved = legacy_mapping[[intensity]]
  } else {
    resolved = .standardizeColnames(intensity)
  }

  if (!(resolved %in% available_cols)) {
    stop(paste0(
      "Intensity column '", intensity, "' not found in input data. ", collapse = ", ")
    )
  }
  resolved
}


#' Assign IsotopeLabelType based on heavy label detection.
#'
#' When \code{heavyLabel} is provided, each row is classified as heavy
#' (\code{"H"}), light (\code{"L"}), or unlabeled (\code{NA}) by inspecting
#' the labeled sequence column for the presence of the label tag.
#'
#' In Spectronaut protein turnover reports, heavy peptides appear in
#' \code{FG.LabeledSequence} with a bracketed modification, e.g.
#' \code{_PEPTIDEK[Lys6]_}.  Any sequence that contains
#' \code{[<heavyLabel>]} is classified as heavy; all others are light.
#' Sequences that do not have amino acids that can carry the label
#' are classified as \code{NA}.  For example, if \code{heavyLabels} is 
#' \code{"Lys6"}, then \code{PEPTIDEZ} is classified as NA since it
#' has no lysine residues that could be labeled.
#' When \code{heavyLabel} is \code{NULL} the column is left untouched so
#' that the downstream \code{columns_to_fill} default of \code{"L"} applies,
#' preserving backwards compatibility.
#'
#' @param spec_input `data.table` after column renaming.
#' @param heavyLabels Character scalar heavy label name (e.g. \code{"Lys6"}),
#'   or \code{NULL}.
#' @return `data.table` with \code{IsotopeLabelType} column added or updated.
#' @keywords internal
#' @noRd
.assignSpectronautIsotopeLabelType = function(spec_input, heavyLabels) {
    PeptideSequence = NULL
    if (is.null(heavyLabels)) {
        return(spec_input)
    }

    bare_amino_acids = sub("\\[.*\\]", "", heavyLabels)
    labeled_aa_regex = paste0(bare_amino_acids, collapse = "|")
    heavy_regex = paste(
        gsub("([\\[\\]])", "\\\\\\1", heavyLabels, perl = TRUE),
        collapse = "|"
    )

    spec_input = .classifyIsotopeLabelType(spec_input, heavy_regex,
                                            labeled_aa_regex = labeled_aa_regex)

    for (i in seq_along(heavyLabels)) {
        escaped = gsub("([\\[\\]])", "\\\\\\1", heavyLabels[i], perl = TRUE)
        spec_input[, PeptideSequence := gsub(escaped, bare_amino_acids[i],
                                             PeptideSequence, perl = TRUE)]
    }
    spec_input
}
