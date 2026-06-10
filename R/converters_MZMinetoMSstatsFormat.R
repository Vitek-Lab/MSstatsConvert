#' Import MZMine files
#'
#' @inheritParams .sharedParametersAmongConverters
#' @inheritParams .cleanRawMZMine
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
#'
#' @details
#' `ProteinName` is assigned from one of three sources, in priority
#' order: the MZMine compound name (mandatory), the SIRIUS name
#' (optional), and an m/z-RT fallback (always available).
#'
#' The **MZMine compound name** is the highest-scoring `compound_name`
#' from `mzmine_annotations` for each feature. This corresponds to MSI
#' Level 2 (Sumner et al. 2007, PMID 27624161): a putative
#' identification by MS/MS spectral matching to a reference library.
#'
#' The **SIRIUS name** comes from SIRIUS's
#' `structure_identifications.tsv` and corresponds to MSI Level 3: an
#' in-silico structure prediction. When `sirius_annotations` is
#' non-NULL, the SIRIUS `name` fills `ProteinName` only for features
#' the MZMine library missed -- the MZMine compound name takes
#' precedence.
#'
#' The **m/z-RT fallback** is an identifier built from the feature's
#' m/z and retention time (for example, `455.282_0.65`). Features that
#' receive no MZMine or SIRIUS annotation are retained, not dropped,
#' and assigned an m/z-RT identifier as their `ProteinName`.
#'
#' Retaining every feature is a deliberate trade-off. A fuller feature
#' set gives more stable medians and a more reliable empirical
#' distribution for global normalization. SIRIUS extends discovery
#' coverage to features that level-2 spectral matching misses. The
#' cost is an increase in the number of hypotheses tested downstream
#' (in `MSstats::groupComparison`), which weakens multiple-testing
#' correction. Users running confirmatory analyses should restrict to
#' the MZMine-annotated features post-conversion; users running
#' discovery analyses benefit from the additional sources despite the
#' FDR burden.
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
#'
#' # With SIRIUS annotations:
#' sirius_path = system.file(
#'   "tinytest/raw_data/MZMine/structure_identifications.tsv",
#'   package = "MSstatsConvert")
#' sirius = data.table::fread(sirius_path)
#' output_with_sirius = MZMinetoMSstatsFormat(
#'   input, annotation = annot,
#'   mzmine_annotations = lib,
#'   sirius_annotations = sirius,
#'   use_log_file = FALSE)
#' head(output_with_sirius)
MZMinetoMSstatsFormat = function(
    input,
    annotation = NULL,
    mzmine_annotations,
    sirius_annotations = NULL,
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

    if (!is.null(sirius_annotations)) {
        sirius_cols = colnames(sirius_annotations)
        missing_sirius = setdiff(c("mappingFeatureId", "name"), sirius_cols)
        if (length(missing_sirius) > 0) {
            stop("sirius_annotations is missing required column(s): ",
                 paste(missing_sirius, collapse = ", "),
                 ". Required: 'mappingFeatureId' and 'name'.")
        }
    }

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "MZMine", ...)
    input = MSstatsConvert::MSstatsClean(
        input,
        mzmine_annotations = mzmine_annotations,
        sirius_annotations = sirius_annotations)
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

    feature_columns = c("PeptideSequence", "PrecursorCharge", "FragmentIon", "ProductCharge")

    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        feature_columns,
        remove_shared_peptides = FALSE,
        remove_single_feature_proteins = FALSE,
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
