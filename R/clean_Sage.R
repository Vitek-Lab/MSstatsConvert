#' Clean raw Sage LFQ output
#'
#' Operates on Sage's `lfq.tsv` report (produced when `quant.lfq: true`). This
#' is a wide-format table with fixed columns `peptide`, `charge`, `proteins`,
#' `q_value`, `score`, `spectral_angle`, followed by one intensity column per
#' input mzML, each headed by the run's file name. Intensity columns are
#' identified as every column that is not one of the six fixed columns (matching
#' is by name, not by file extension, so renamed files are handled). The table is
#' melted to long format, columns are renamed to the MSstats standard, and zero
#' intensities -- which Sage writes for precursors it did not quantify in a run --
#' are converted to `NA`.
#'
#' @param msstats_object an object of class `MSstatsSageFiles`.
#' @return data.table
#' @keywords internal
.cleanRawSage = function(msstats_object) {
    Intensity = NULL

    sage_input = getInputFile(msstats_object, "input")
    sage_input = data.table::as.data.table(sage_input)

    fixed_columns = c("peptide", "charge", "proteins", "q_value",
                      "score", "spectral_angle")
    required_columns = c("peptide", "charge", "proteins", "q_value")
    missing_columns = setdiff(required_columns, colnames(sage_input))
    if (length(missing_columns) > 0) {
        msg = paste("The following required columns are missing from the Sage",
                    "input:", paste(missing_columns, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }

    intensity_columns = setdiff(colnames(sage_input), fixed_columns)
    if (length(intensity_columns) == 0) {
        msg = paste("No intensity columns found in the Sage input. Expected at",
                    "least one per-run intensity column in addition to the fixed",
                    "columns:", paste(fixed_columns, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }

    id_columns = intersect(c("proteins", "peptide", "charge", "q_value"),
                           colnames(sage_input))
    sage_input = sage_input[, c(id_columns, intensity_columns), with = FALSE]

    long = data.table::melt(sage_input,
                            id.vars = id_columns,
                            measure.vars = intensity_columns,
                            variable.name = "Run",
                            value.name = "Intensity",
                            variable.factor = FALSE)

    data.table::setnames(long,
                         c("proteins", "peptide", "charge"),
                         c("ProteinName", "PeptideSequence", "PrecursorCharge"))

    long[, Intensity := as.numeric(Intensity)]
    long[Intensity == 0, Intensity := NA_real_]

    if (all(long$PrecursorCharge == -1)) {
        msg = paste("** All PrecursorCharge values are -1: Sage combined charge",
                    "states (combine_charge_states = true), so the feature key is",
                    "effectively the peptide sequence alone.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
    }

    .logSuccess("Sage", "clean")
    long
}
