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
    PrecursorCharge = FragmentIon = ProductCharge = IsotopeLabelType = NULL
    sample_col = id = score = compound_name = NULL

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

    if (!is.null(mzmine_annotations)) {
        ann <- data.table::as.data.table(mzmine_annotations)
        required_ann <- c("id", "compound_name", "score")
        missing_ann <- setdiff(required_ann, colnames(ann))
        if (length(missing_ann) > 0) {
            stop("mzmine_annotations is missing required column(s): ",
                 paste(missing_ann, collapse = ", "), ".")
        }
        data.table::setorder(ann, id, -score)
        ann_top <- unique(ann, by = "id")
        matched <- ann_top[match(mz_input[[id_col]], ann_top[["id"]]),
                           compound_name]
        compound <- ifelse(is.na(matched), mz_rt_fallback, matched)
    } else {
        compound <- mz_rt_fallback
    }

    mz_input[, ProteinName := compound]
    mz_input[, PeptideSequence := as.character(get(id_col))]

    long <- data.table::melt(
        mz_input,
        id.vars = c("ProteinName", "PeptideSequence"),
        measure.vars = peak_area_cols,
        variable.name = "sample_col",
        value.name = "Intensity",
        variable.factor = FALSE)

    long[, PrecursorCharge := NA_integer_]
    long[, FragmentIon := NA_character_]
    long[, ProductCharge := NA_integer_]
    long[, IsotopeLabelType := "Light"]
    long[, Run := sub(paste0(peak_area_suffix, "$"), "", sample_col)]

    final_cols <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
                    "FragmentIon", "ProductCharge", "IsotopeLabelType",
                    "Run", "Intensity")
    long <- long[, final_cols, with = FALSE]

    .logSuccess("MZMine", "clean")
    long
}
