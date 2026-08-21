# Test .cleanRawSage on the Sage lfq.tsv fixture ---------------------------
sage_path = system.file("tinytest/raw_data/Sage/lfq.tsv",
                        package = "MSstatsConvert")
if (!nzchar(sage_path)) {
    exit_file("Sage fixtures not present in inst/tinytest/raw_data/Sage/")
}
sage_raw = data.table::fread(sage_path)

fixed_cols = c("peptide", "charge", "proteins", "q_value", "score",
               "spectral_angle")
intensity_cols = setdiff(colnames(sage_raw), fixed_cols)

msstats_input = MSstatsConvert::MSstatsImport(
    list(input = sage_raw), "MSstats", "Sage")
cleaned = MSstatsConvert:::.cleanRawSage(msstats_input)
cleaned = data.table::as.data.table(cleaned)

# Renamed to canonical MSstats names; wide reshaped to long
expect_true(all(c("ProteinName", "PeptideSequence", "PrecursorCharge",
                  "q_value", "Run", "Intensity") %in% colnames(cleaned)))
expect_false("proteins" %in% colnames(cleaned))
expect_false("peptide" %in% colnames(cleaned))
expect_false("charge" %in% colnames(cleaned))
# score / spectral_angle are dropped in cleaning
expect_false("score" %in% colnames(cleaned))
expect_false("spectral_angle" %in% colnames(cleaned))

# One Run per intensity column; every wide row melted across every run
expect_equal(data.table::uniqueN(cleaned$Run), length(intensity_cols))
expect_equal(nrow(cleaned), nrow(sage_raw) * length(intensity_cols))

# Zero intensities converted to NA; NA count equals the number of zero cells
expect_false(any(cleaned$Intensity == 0, na.rm = TRUE))
n_zero_cells = sum(as.matrix(sage_raw[, intensity_cols, with = FALSE]) == 0,
                   na.rm = TRUE)
expect_equal(sum(is.na(cleaned$Intensity)), n_zero_cells)

# q_value is retained for the downstream converter-level filter
expect_true("q_value" %in% colnames(cleaned))

# Combined charge states -> PrecursorCharge is -1 throughout on this fixture
expect_true(all(cleaned$PrecursorCharge == -1))

# A missing required column is an error naming the offending column
bad_input = data.table::copy(sage_raw)
bad_input$q_value = NULL
msstats_bad = MSstatsConvert::MSstatsImport(
    list(input = bad_input), "MSstats", "Sage")
expect_error(MSstatsConvert:::.cleanRawSage(msstats_bad), "q_value")

# No intensity columns is an error
only_fixed = sage_raw[, intersect(fixed_cols, colnames(sage_raw)), with = FALSE]
msstats_only_fixed = MSstatsConvert::MSstatsImport(
    list(input = only_fixed), "MSstats", "Sage")
expect_error(MSstatsConvert:::.cleanRawSage(msstats_only_fixed), "intensity")


# Test .cleanRawSage on the charge-resolved fixture ------------------------
sage_cr_path = system.file("tinytest/raw_data/Sage/lfq_charge_resolved.tsv",
                           package = "MSstatsConvert")
if (nzchar(sage_cr_path)) {
    sage_cr = data.table::fread(sage_cr_path)
    msstats_cr = MSstatsConvert::MSstatsImport(
        list(input = sage_cr), "MSstats", "Sage")
    cleaned_cr = MSstatsConvert:::.cleanRawSage(msstats_cr)
    cleaned_cr = data.table::as.data.table(cleaned_cr)

    # Real precursor charges survive cleaning (not all -1)
    expect_true(all(c(2, 3) %in% unique(cleaned_cr$PrecursorCharge)))
    expect_false(all(cleaned_cr$PrecursorCharge == -1))
}
