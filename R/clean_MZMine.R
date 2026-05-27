#' Clean raw MZMine files
#'
#' Operates on the column names produced by MZMine after MSstatsConvert's
#' internal column-name standardization (spaces collapsed and dots removed):
#' "row ID" becomes `rowID`, and each "<sample> Peak area" becomes
#' `<standardized-sample>Peakarea`.
#'
#' @param msstats_object an object of class `MSstatsMZMineFiles`.
#' @param mzmine_annotations `data.frame` of MZMine spectral-library
#'   annotations with columns `id`, `compound_name`, `score`. Required;
#'   passing `NULL` raises an error. The highest-scoring `compound_name`
#'   per feature is used as `ProteinName`, and features in the quant
#'   table with no matching annotation row are dropped from the output.
#'   These are MSI Level 2 annotations (putative identification via
#'   MS/MS spectral matching). See the public `MZMinetoMSstatsFormat`
#'   docstring for the full scope discussion.
#' @return data.table
#' @keywords internal
.cleanRawMZMine <- function(msstats_object, mzmine_annotations) {
    ProteinName = PeptideSequence = Intensity = Run = NULL
    PrecursorCharge = FragmentIon = ProductCharge = NULL
    id = score = compound_name = i.compound_name = NULL

    mz_input = getInputFile(msstats_object, "input")
    mz_input = data.table::as.data.table(mz_input)

    peak_area_suffix <- "Peakarea"
    peak_area_cols <- grep(paste0(peak_area_suffix, "$"),
                           colnames(mz_input), value = TRUE)
    if (length(peak_area_cols) == 0) {
        stop("No 'Peak area' columns found in the input. Expected per-sample ",
             "columns named '<run> Peak area' (e.g. 'sampleA.mzML Peak area').")
    }
    id_col <- "rowID"
    required_meta <- id_col
    missing_meta <- setdiff(required_meta, colnames(mz_input))
    if (length(missing_meta) > 0) {
        stop("Missing required MZMine metadata column (expected 'row ID'). ",
             "After standardization, looked for: ",
             paste(missing_meta, collapse = ", "), ".")
    }

    if (is.null(mzmine_annotations)) {
        stop("mzmine_annotations is required. Pass a data.frame with ",
             "columns 'id', 'compound_name', 'score'.")
    }
    feature_to_compound <- data.table::as.data.table(mzmine_annotations)
    required_ann <- c("id", "compound_name", "score")
    missing_ann <- setdiff(required_ann, colnames(feature_to_compound))
    if (length(missing_ann) > 0) {
        stop("mzmine_annotations is missing required column(s): ",
             paste(missing_ann, collapse = ", "), ".")
    }
    feature_to_compound[, score := suppressWarnings(as.numeric(score))]
    if (anyNA(feature_to_compound$score)) {
        stop("The 'score' column in the mzmine annotations file must contain numeric values.")
    }
    data.table::setorder(feature_to_compound, id, -score)
    feature_to_compound <- unique(feature_to_compound, by = "id")
    # Inner-join filter: drop quant rows with no matching annotation.
    mz_input[
        feature_to_compound,
        ProteinName := i.compound_name,
        on = setNames("id", id_col)
    ]
    mz_input <- mz_input[!is.na(ProteinName)]

    retained_ids <- feature_to_compound$id
    retained_msg <- paste0("** MZMine: retained ", length(retained_ids),
                           " feature(s) after annotation join: ",
                           paste(retained_ids, collapse = ", "))
    getOption("MSstatsLog")("INFO", retained_msg)
    getOption("MSstatsMsg")("INFO", retained_msg)

    mz_input[, PeptideSequence := as.character(get(id_col))]

    long <- data.table::melt(
        mz_input,
        id.vars = c("ProteinName", "PeptideSequence"),
        measure.vars = peak_area_cols,
        variable.name = "Run",
        value.name = "Intensity",
        variable.factor = FALSE)

    long[, PrecursorCharge := NA_integer_]
    long[, FragmentIon := NA_character_]
    long[, ProductCharge := NA_integer_]
    long[, Run := sub(paste0(peak_area_suffix, "$"), "", Run)]

    final_cols <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
                    "FragmentIon", "ProductCharge",
                    "Run", "Intensity")
    long <- long[, final_cols, with = FALSE]

    .logSuccess("MZMine", "clean")
    long
}
