#' Clean raw MZMine files
#'
#' Operates on the column names produced by MZMine after MSstatsConvert's
#' internal column-name standardization (spaces collapsed and dots removed):
#' "row ID" becomes `rowID`, "row m/z" becomes `rowmz`, "row retention time"
#' becomes `rowretentiontime`, and each "<sample> Peak area" becomes
#' `<standardized-sample>Peakarea`.
#'
#' @param msstats_object an object of class `MSstatsMZMineFiles`.
#' @param mzmine_annotations optional `data.frame` of MZMine spectral-library
#'   annotations with columns `id`, `compound_name`, `score`. When supplied,
#'   the highest-scoring `compound_name` per feature is used as `ProteinName`.
#'   Features without a matching annotation row fall back to an mz_rt string
#'   `paste0(round(mz, 4), "_", round(rt, 2))`. When `NULL`, every feature
#'   uses the mz_rt fallback.
#' @return data.table
#' @keywords internal
.cleanRawMZMine <- function(msstats_object, mzmine_annotations = NULL) {
    ProteinName = PeptideSequence = Intensity = Run = NULL
    PrecursorCharge = FragmentIon = ProductCharge = NULL
    id = score = compound_name = i.compound_name = NULL

    mz_input <- getInputFile(msstats_object, "input")
    mz_input <- data.table::as.data.table(mz_input)

    peak_area_suffix <- "Peakarea"
    peak_area_cols <- grep(paste0(peak_area_suffix, "$"),
                           colnames(mz_input), value = TRUE)
    if (length(peak_area_cols) == 0) {
        stop("No 'Peak area' columns found in the input. Expected per-sample ",
             "columns named '<run> Peak area' (e.g. 'sampleA.mzML Peak area').")
    }
    id_col <- "rowID"
    mz_col <- "rowmz"
    rt_col <- "rowretentiontime"
    required_meta <- c(id_col, mz_col, rt_col)
    missing_meta <- setdiff(required_meta, colnames(mz_input))
    if (length(missing_meta) > 0) {
        stop("Missing required MZMine metadata column(s) (expected 'row ID', ",
             "'row m/z', 'row retention time'). After standardization, ",
             "looked for: ", paste(missing_meta, collapse = ", "), ".")
    }

    mz_rt_fallback <- paste0(round(mz_input[[mz_col]], 4), "_",
                             round(mz_input[[rt_col]], 2))
    mz_input[, ProteinName := mz_rt_fallback]

    if (!is.null(mzmine_annotations)) {
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
        # Sort by id ascending and score descending so the highest-scoring
        # annotation per id is the first row in each group.
        data.table::setorder(feature_to_compound, id, -score)
        # Collapse to one row per id (the highest-scoring). data.table's
        # unique() with a 'by' arg keeps the first row per group, which after
        # the sort above is the highest-scoring annotation.
        feature_to_compound <- unique(feature_to_compound, by = "id")
        # Join: unmatched mz_input rows keep the mz_rt_fallback ProteinName
        # set above.
        mz_input[
            feature_to_compound,
            ProteinName := i.compound_name,
            on = setNames("id", id_col)
        ]
    }

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
