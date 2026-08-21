#' Import Sage LFQ files
#'
#' Converts the label-free quantification report (`lfq.tsv`) produced by the Sage
#' search engine into a `data.frame` in the format required by MSstats. The input
#' is wide (one row per precursor, one intensity column per run); it is reshaped
#' to long format, q-value filtered, and returned ready for `dataProcess`.
#'
#' @inheritParams .sharedParametersAmongConverters
#' @param input Sage `lfq.tsv` report, as a `data.frame`/`data.table` or a path.
#'   Wide format with fixed columns `peptide`, `charge`, `proteins`, `q_value`,
#'   `score`, `spectral_angle`, followed by one intensity column per input mzML,
#'   each headed by the run's file name.
#' @param annotation `data.frame` with `Run`, `Condition` and `BioReplicate`
#'   columns (a `Fraction` column may also be supplied). This argument is
#'   **required**: Sage's `lfq.tsv` carries no experimental design, so condition
#'   and replicate information must be provided separately.
#' @param qvalue_cutoff Cutoff for the `q_value` column. Default is 0.01.
#' @param filter_with_Qvalue TRUE (default) replaces intensities whose `q_value`
#'   is greater than `qvalue_cutoff` with `NA` (treated as censored missing
#'   downstream). FALSE performs no q-value filtering. Sage does not pre-filter
#'   `lfq.tsv` on `lfq_settings.peptide_q_value`, so this filter is load-bearing.
#'
#' @return `data.frame` in the MSstats required format.
#'
#' @section Input file:
#' Use `lfq.tsv`, which Sage writes only when `quant.lfq: true` is set in the
#' search configuration. Do **not** use `results.sage.tsv`: its `ms2_intensity`
#' column is the summed intensity of matched b/y fragment ions -- a PSM score
#' feature -- and is not a quantitative measure of precursor abundance.
#'
#' @section Charge states:
#' Sage's `combine_charge_states` option (default `true`) sums charge states and
#' writes `charge` as `-1` for every row, so `PrecursorCharge` will be `-1`
#' throughout and the feature key reduces to the peptide sequence. Setting
#' `combine_charge_states: false` reports real precursor charges, but is
#' considerably slower across multiple files.
#'
#' @section Run name matching:
#' MSstatsConvert standardizes column names by removing spaces and dots (`.`)
#' while preserving hyphens and underscores. The melted `Run` values (the
#' intensity column headers) and the annotation `Run` values are both passed
#' through this same standardization before merging, so they match automatically.
#' For example, a run named `B.naive_01steady-state.mzML.gz` in the annotation
#' becomes `Bnaive_01steady-statemzMLgz`; supply the raw file name in the
#' annotation and the merge resolves it. Note that the `Run` values in the
#' returned table are the standardized form.
#'
#' @section Shared peptides:
#' Sage pre-joins shared proteins into a single semicolon-delimited `proteins`
#' value (e.g. `sp|A|X;sp|B|Y`). Because that is one `ProteinName` string rather
#' than several, MSstats' shared-peptide removal sees a single protein per
#' peptide and does not treat these rows as shared. Consequently
#' `useUniquePeptide = TRUE` has no effect on peptides that Sage reports against a
#' shared (semicolon-joined) protein group.
#'
#' @export
#'
#' @examples
#' sage_lfq = system.file("tinytest/raw_data/Sage/lfq.tsv",
#'                        package = "MSstatsConvert")
#' annot_path = system.file("tinytest/raw_data/Sage/annotation.csv",
#'                          package = "MSstatsConvert")
#' if (nzchar(sage_lfq) && nzchar(annot_path)) {
#'     sage_input = data.table::fread(sage_lfq)
#'     annotation = read.csv(annot_path)
#'     sage_imported = SagetoMSstatsFormat(sage_input, annotation,
#'                                         use_log_file = FALSE)
#'     head(sage_imported)
#' }
#'
SagetoMSstatsFormat = function(
        input, annotation, useUniquePeptide = TRUE,
        removeFewMeasurements = TRUE, removeProtein_with1Peptide = FALSE,
        qvalue_cutoff = 0.01, filter_with_Qvalue = TRUE,
        use_log_file = TRUE, append = FALSE, verbose = TRUE, log_file_path = NULL,
        ...
) {
    IsotopeLabelType = NULL

    validation_config = list(
        input = input,
        annotation = annotation,
        filter_with_Qvalue = filter_with_Qvalue,
        qvalue_cutoff = qvalue_cutoff,
        useUniquePeptide = useUniquePeptide,
        removeFewMeasurements = removeFewMeasurements,
        removeProtein_with1Feature = removeProtein_with1Peptide,
        use_log_file = use_log_file,
        append = append,
        verbose = verbose,
        log_file_path = log_file_path
    )
    .validateMSstatsConverterParameters(validation_config)

    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstats", "Sage", ...)
    input = MSstatsConvert::MSstatsClean(input)

    if (inherits(annotation, "data.frame") &&
        is.element("IsotopeLabelType", colnames(annotation))) {
        annotation = data.table::as.data.table(annotation)
        annotation[, IsotopeLabelType := NULL]
        msg = paste("** An IsotopeLabelType column was found in the annotation",
                    "and has been dropped. Sage LFQ is label-free;",
                    "IsotopeLabelType is set to 'L' for all rows.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, annotation)

    qval_filter = list(score_column = "q_value",
                       score_threshold = qvalue_cutoff,
                       direction = "smaller",
                       behavior = "fill",
                       handle_na = "keep",
                       fill_value = NA_real_,
                       filter = filter_with_Qvalue,
                       drop_column = TRUE)

    feature_columns = c("PeptideSequence", "PrecursorCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        feature_columns,
        remove_shared_peptides = useUniquePeptide,
        remove_single_feature_proteins = removeProtein_with1Peptide,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = function(x, na.rm = TRUE) {
                if (all(is.na(x))) NA_real_ else max(x, na.rm = na.rm)
            }),
        score_filtering = list(qvalue = qval_filter),
        columns_to_fill = list("FragmentIon" = NA,
                               "ProductCharge" = NA,
                               "IsotopeLabelType" = "L"))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns,
                                                  remove_few = removeFewMeasurements)

    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}
