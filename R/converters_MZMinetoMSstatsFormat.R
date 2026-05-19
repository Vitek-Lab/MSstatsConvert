#' Import MZMine files
#'
#' @inheritParams .sharedParametersAmongConverters
#' @param input MZMine feature-quantification table (wide format; one row per
#'   feature). Must include the metadata columns `row ID`, `row m/z`,
#'   `row retention time`, and per-sample peak-area columns named
#'   `"<run> Peak area"` (e.g. `"sampleA.mzML Peak area"`).
#' @param annotation `data.frame` with columns `Run`, `Condition`,
#'   `BioReplicate`. `Run` values must match the sample column names with the
#'   trailing `" Peak area"` stripped.
#' @param mzmine_annotations optional `data.frame` of MZMine spectral-library
#'   annotations with columns `id`, `compound_name`, `score`. When supplied,
#'   the highest-scoring `compound_name` per feature is used as `ProteinName`;
#'   features without a matching annotation row fall back to an mz_rt string
#'   `paste0(round(mz, 4), "_", round(rt, 2))`. When `NULL`, every feature
#'   uses the mz_rt fallback.
#' @param removeProtein_with1Feature `TRUE` will remove proteins (compounds)
#'   represented by a single feature. Default `FALSE`.
#' @param summaryforMultipleRows `max` (default) or `sum` — used when multiple
#'   rows map to the same feature/run combination.
#'
#' @return data.table in the MSstats required format.
#'
#' @export
#'
#' @examples
#' input_path = system.file("tinytest/raw_data/MZMine/mzmine_input.csv",
#'                          package = "MSstatsConvert")
#' annot_path = system.file("tinytest/raw_data/MZMine/annotation.csv",
#'                          package = "MSstatsConvert")
#' lib_path   = system.file("tinytest/raw_data/MZMine/mzmine_annotations.csv",
#'                          package = "MSstatsConvert")
#' input = data.table::fread(input_path)
#' annot = data.table::fread(annot_path)
#' lib   = data.table::fread(lib_path)
#' output = MZMinetoMSstatsFormat(input, annotation = annot,
#'                                mzmine_annotations = lib,
#'                                use_log_file = FALSE)
#' head(output)
MZMinetoMSstatsFormat = function(
    input,
    annotation = NULL,
    mzmine_annotations = NULL,
    removeProtein_with1Feature = FALSE,
    summaryforMultipleRows = max,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL,
    ...) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "MZMine")
    input = MSstatsConvert::MSstatsClean(
        input, mzmine_annotations = mzmine_annotations)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
    fill_isotope_label_type = if ("IsotopeLabelType" %in% colnames(input))
        list() else list("IsotopeLabelType" = "Light")

    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        feature_columns,
        remove_shared_peptides = FALSE,
        remove_single_feature_proteins = removeProtein_with1Feature,
        exact_filtering = NULL,
        pattern_filtering = NULL,
        aggregate_isotopic = FALSE,
        feature_cleaning = list(
            remove_features_with_few_measurements = FALSE,
            summarize_multiple_psms = summaryforMultipleRows),
        columns_to_fill = c(list(Fraction = 1), fill_isotope_label_type))
    input[, Intensity := ifelse(Intensity == 0, NA, Intensity)]

    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  fill_incomplete = TRUE,
                                                  handle_fractions = FALSE,
                                                  remove_few = FALSE)

    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}
