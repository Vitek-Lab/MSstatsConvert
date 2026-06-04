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
#'   per feature is used as `ProteinName` (tier 1, MSI Level 2 putative
#'   identification via MS/MS spectral matching). See the public
#'   `MZMinetoMSstatsFormat` docstring for the full tier discussion.
#' @param sirius_annotations Optional `data.frame` of SIRIUS
#'   `structure_identifications.tsv` output, or `NULL`. Only the
#'   `mappingFeatureId` and `name` columns are read. When supplied,
#'   SIRIUS `name` fills `ProteinName` for features that received no
#'   MZMine compound (tier 2, MSI Level 3).
#' @return data.table
#' @keywords internal
.cleanRawMZMine <- function(msstats_object, mzmine_annotations,
                            sirius_annotations = NULL) {
    ProteinName = PeptideSequence = Intensity = Run = NULL
    PrecursorCharge = FragmentIon = ProductCharge = NULL
    id = score = compound_name = i.compound_name = NULL
    rowmz = rowretentiontime = mappingFeatureId = name = NULL

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
    mz_col <- "rowmz"
    rt_col <- "rowretentiontime"
    required_meta <- c(id_col, mz_col, rt_col)
    missing_meta <- setdiff(required_meta, colnames(mz_input))
    if (length(missing_meta) > 0) {
        stop("Missing required MZMine metadata column(s) ",
             "(expected 'row ID', 'row m/z', 'row retention time'). ",
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
    # Tier 1: MZMine compound name (left-join, no drop).
    mz_input[
        feature_to_compound,
        ProteinName := i.compound_name,
        on = setNames("id", id_col)
    ]
    n_tier1 <- sum(!is.na(mz_input$ProteinName))

    # Tier 2: SIRIUS name fills features still NA after tier 1.
    n_tier2 <- 0L
    if (!is.null(sirius_annotations)) {
        sirius_dt <- data.table::as.data.table(sirius_annotations)
        sirius_dt <- sirius_dt[, c("mappingFeatureId", "name"), with = FALSE]
        sirius_dt[, name := ifelse(is.na(name) | name == "",
                                   NA_character_, as.character(name))]
        sirius_dt[, mappingFeatureId := as.character(mappingFeatureId)]
        data.table::setorder(sirius_dt, mappingFeatureId)
        sirius_dt <- unique(sirius_dt, by = "mappingFeatureId")
        mz_input[, ProteinName := ifelse(
            is.na(ProteinName),
            sirius_dt[
                .(as.character(get(id_col))),
                on = "mappingFeatureId",
                name
            ],
            ProteinName)]
        n_tier2 <- sum(!is.na(mz_input$ProteinName)) - n_tier1
    }

    # Tier 3: mz_rt fallback for features still NA.
    na_mask <- is.na(mz_input$ProteinName)
    n_tier3 <- sum(na_mask)
    if (n_tier3 > 0) {
        mz_input[na_mask, ProteinName := paste0(
            round(get(mz_col), 4), "_", round(get(rt_col), 2))]
    }

    tier_msg <- paste0(
        "** MZMine ProteinName assignment: ",
        "tier 1 (MZMine compound): ", n_tier1, " feature(s); ",
        "tier 2 (SIRIUS name): ", n_tier2, " feature(s); ",
        "tier 3 (mz_rt fallback): ", n_tier3, " feature(s).")
    getOption("MSstatsLog")("INFO", tier_msg)
    getOption("MSstatsMsg")("INFO", tier_msg)

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
