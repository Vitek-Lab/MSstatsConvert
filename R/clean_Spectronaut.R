#' Clean raw Spectronaut output.
#' @param msstats_object an object of class `MSstatsSpectronautFiles`.
#' @inheritParams SpectronauttoMSstatsFormat
#' @return `data.table`
#' @keywords internal
.cleanRawSpectronaut = function(msstats_object, intensity,
                                calculateAnomalyScores,
                                anomalyModelFeatures,
                                heavyLabel = NULL,
                                labelColumn = "FG.LabeledSequence") {
  FFrgLossType = FExcludedFromQuantification = NULL

  spec_input = getInputFile(msstats_object, "input")
  .validateSpectronautInput(spec_input)

  # --- Normalize missing columns that vary across Spectronaut report formats ---
  spec_input = .addMissingSpectronautColumns(spec_input)

  spec_input = spec_input[FFrgLossType == "noloss", ]

  f_charge_col = .findAvailable(c("FCharge", "FFrgZ"), colnames(spec_input))
  pg_qval_col = .findAvailable(c("PGQvalue"), colnames(spec_input))
  interference_col = .findAvailable(c("FPossibleInterference"),
                                    colnames(spec_input))
  exclude_col = .findAvailable(c("FExcludedFromQuantification"),
                               colnames(spec_input))

  # Resolve intensity column: accepts enum alias OR raw standardized column name
  intensity_column = .resolveSpectronautIntensityColumn(intensity, colnames(spec_input))

  # Resolve peptide sequence column: prefer EGModifiedSequence, fall back to
  # FGLabeledSequence (protein turnover format uses FG.LabeledSequence)
  peptide_col = .findAvailable(c("EGModifiedSequence", "FGLabeledSequence"),
                               colnames(spec_input))

  # Resolve protein name column: prefer PGProteinGroups, fall back to
  # PGProteinAccessions (protein turnover format omits PG.ProteinGroups)
  protein_col = .findAvailable(c("PGProteinGroups", "PGProteinAccessions"),
                               colnames(spec_input))

  # Resolve BioReplicate column: prefer RReplicate, fall back to RCondition
  # (protein turnover format does not include R.Replicate)
  replicate_col = .findAvailable(c("RReplicate", "RCondition"),
                                 colnames(spec_input))

  cols = c(protein_col, peptide_col, "FGCharge", "FFrgIon",
           f_charge_col, "RFileName", "RCondition", replicate_col,
           "EGQvalue", pg_qval_col, interference_col, exclude_col,
           intensity_column)
  if (calculateAnomalyScores){
    cols = c(cols, anomalyModelFeatures)
  }
  cols = intersect(cols, colnames(spec_input))
  spec_input = spec_input[, cols, with = FALSE]

  data.table::setnames(
    spec_input,
    c(protein_col, peptide_col, "FGCharge", "FFrgIon",
      f_charge_col, "RFileName", intensity_column,
      "RCondition", replicate_col),
    c("ProteinName", "PeptideSequence", "PrecursorCharge", "FragmentIon",
      "ProductCharge", "Run", "Intensity", "Condition", "BioReplicate"),
    skip_absent = TRUE)

  # Assign IsotopeLabelType based on heavy label detection when requested
  spec_input = .assignSpectronautIsotopeLabelType(
    spec_input, heavyLabel, labelColumn, msstats_object)

  .logSuccess("Spectronaut", "clean")
  spec_input
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
.addMissingSpectronautColumns = function(spec_input) {
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
.resolveSpectronautIntensityColumn = function(intensity, available_cols) {
  legacy_mapping = c(
    "PeakArea"           = "FPeakArea",
    "NormalizedPeakArea" = "FNormalizedPeakArea",
    "MS1Quantity"        = "FGMS1Quantity",
    "MS2Quantity"        = "FGMS2Quantity"
  )

  if (intensity %in% names(legacy_mapping)) {
    resolved = legacy_mapping[[intensity]]
  } else {
    # Treat as a raw standardized column name
    resolved = intensity
  }

  if (!(resolved %in% available_cols)) {
    stop(paste0(
      "Intensity column '", resolved, "' not found in input data. ",
      "Available columns include: ",
      paste(grep("Quantity|PeakArea", available_cols, value = TRUE), collapse = ", ")
    ))
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
#' Sequences that belong to peptide families that cannot carry the label
#' (i.e. the same stripped sequence never appears in a heavy form in the
#' entire dataset) are classified as \code{NA}.
#'
#' When \code{heavyLabel} is \code{NULL} the column is left untouched so
#' that the downstream \code{columns_to_fill} default of \code{"L"} applies,
#' preserving backwards compatibility.
#'
#' @param spec_input `data.table` after column renaming.
#' @param heavyLabel Character scalar heavy label name (e.g. \code{"Lys6"}),
#'   or \code{NULL}.
#' @param labelColumn Raw (dot-separated) column name that holds the labeled
#'   sequence (e.g. \code{"FG.LabeledSequence"}).
#' @param msstats_object The original MSstats object (used to access the
#'   standardized label column after import).
#' @return `data.table` with \code{IsotopeLabelType} column added or updated.
#' @keywords internal
.assignSpectronautIsotopeLabelType = function(spec_input, heavyLabel,
                                              labelColumn, msstats_object) {
  IsotopeLabelType = PeptideSequence = NULL

  if (is.null(heavyLabel)) {
    return(spec_input)
  }

  # The label column may have already been renamed to PeptideSequence if it was
  # the chosen peptide column.  We need the original labeled sequence values.
  # Retrieve them from the cleaned input (PeptideSequence column).
  if (!("PeptideSequence" %in% colnames(spec_input))) {
    msg = paste0("Cannot assign IsotopeLabelType: 'PeptideSequence' column ",
                 "not found after cleaning. Skipping label assignment.")
    getOption("MSstatsLog")("WARN", msg)
    getOption("MSstatsMsg")("WARN", msg)
    return(spec_input)
  }

  heavy_pattern = paste0("[", heavyLabel, "]")

  spec_input[, IsotopeLabelType := data.table::fifelse(
    grepl(heavy_pattern, PeptideSequence, fixed = TRUE),
    "H",
    "L"
  )]

  # Identify stripped sequences that appear ONLY as light (no heavy counterpart
  # exists anywhere in the dataset).  These peptides cannot be labelled and
  # should receive NA rather than "L" to distinguish them from the light
  # channel of a quantified heavy/light pair.
  stripped_col = .findAvailable(
    c("PEPStrippedSequence", "PeptideSequence"), colnames(spec_input))

  if (!is.null(stripped_col) && stripped_col != "PeptideSequence") {
    heavy_sequences = spec_input[IsotopeLabelType == "H",
                                 unique(get(stripped_col))]
    spec_input[IsotopeLabelType == "L" &
                 !(get(stripped_col) %in% heavy_sequences),
               IsotopeLabelType := NA_character_]
  }

  msg = paste0("** IsotopeLabelType assigned using heavy label: '", heavyLabel,
               "'. Heavy (H): ",
               sum(spec_input$IsotopeLabelType == "H", na.rm = TRUE),
               ", Light (L): ",
               sum(spec_input$IsotopeLabelType == "L", na.rm = TRUE),
               ", Unlabeled (NA): ",
               sum(is.na(spec_input$IsotopeLabelType)))
  getOption("MSstatsLog")("INFO", msg)
  getOption("MSstatsMsg")("INFO", msg)

  spec_input
}


#' Helper method to validate input has necessary columns
#' @param spec_input dataframe input
#' @noRd
.validateSpectronautInput = function(spec_input) {
    # Only FGCharge is truly required; all other formerly-required columns are
    # either synthesized by .addMissingSpectronautColumns or detected via
    # .findAvailable fallbacks so that protein turnover reports (which omit
    # several standard columns) are handled without pre-processing by the caller.
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
}
