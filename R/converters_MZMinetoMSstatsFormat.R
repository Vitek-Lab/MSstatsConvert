#' Import MZMine files
#'
#' @inheritParams .sharedParametersAmongConverters
#' @param input MZMine feature-quantification table (wide format; one row per
#'   feature). Must include the metadata columns `row ID`, `row m/z`,
#'   `row retention time`, and per-sample peak-area columns named
#'   `"<run> Peak area"` (e.g. `"sampleA.mzML Peak area"`).
#' @param annotation `data.frame` with columns `Run`, `Condition`,
#'   `BioReplicate`. `Run` values must match MSstatsConvert-standardized sample
#'   names (after column-name normalization removes spaces and dots) with the
#'   trailing `"Peakarea"` suffix removed. For example, a quant-file column
#'   `"sampleA.mzML Peak area"` becomes `"sampleAmzML"` after standardization,
#'   so the corresponding `Run` value must be `sampleAmzML`.
#' @param mzmine_annotations `data.frame` of MZMine spectral-library
#'   annotations with columns `id`, `compound_name`, `score`. Required:
#'   the highest-scoring `compound_name` per feature is used as
#'   `ProteinName`, and features in the quant table with no matching
#'   annotation row are dropped from the output.
#'
#'   These are MSI Level 2 annotations (putative identification via
#'   MS/MS spectral matching against a reference library). Higher-
#'   confidence Level 1 identifications require pure reference standards
#'   and are out of scope here. Lower-confidence annotations such as
#'   Level 3 (SIRIUS, MS2Query) or Level 4 (molecular formula via
#'   CANOPUS) are not currently supported -- features without a Level 2
#'   annotation row are filtered out.
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
    mzmine_annotations,
    removeProtein_with1Feature = FALSE,
    summaryforMultipleRows = max,
    use_log_file = TRUE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL,
    ...) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)

    if (missing(mzmine_annotations) || is.null(mzmine_annotations)) {
        stop("mzmine_annotations is required. Pass a data.frame with ",
             "columns 'id', 'compound_name', 'score'.")
    }

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "MZMine", ...)
    input = MSstatsConvert::MSstatsClean(
        input, mzmine_annotations = mzmine_annotations)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")

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
        columns_to_fill = list(Fraction = 1, IsotopeLabelType = "Light"))
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
