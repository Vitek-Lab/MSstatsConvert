# Tests for SpectronauttoMSstatsFormat with protein turnover (BoxCar) data.
#
# The BoxCar report differs from a standard Spectronaut export in several ways:
#   - No PG.ProteinGroups  (uses PG.ProteinAccessions instead)
#   - No F.FrgLossType     (synthesized as "noloss")
#   - No F.ExcludedFromQuantification (synthesized as FALSE)
#   - No R.Replicate       (falls back to R.Condition for BioReplicate)
#   - No EG.ModifiedSequence (uses FG.LabeledSequence)
#   - No F.FrgIon / F.Charge (synthesized as NA)
#   - Intensity sourced from FG.MS1Quantity or FG.MS2Quantity
#   - Heavy peptides identified by a bracketed label, e.g. [Lys6]

boxcar_path = system.file(
    "tinytest/raw_data/Spectronaut/boxcar_protein_turnover_input.csv",
    package = "MSstatsConvert")
boxcar_raw = data.table::fread(boxcar_path)


# --- Basic format conversion (no heavy label) --------------------------------

output_basic = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "MS1Quantity",
    use_log_file = FALSE
)

expect_true("Run"              %in% colnames(output_basic))
expect_true("ProteinName"      %in% colnames(output_basic))
expect_true("PeptideSequence"  %in% colnames(output_basic))
expect_true("PrecursorCharge"  %in% colnames(output_basic))
expect_true("Intensity"        %in% colnames(output_basic))
expect_true("FragmentIon"      %in% colnames(output_basic))
expect_true("ProductCharge"    %in% colnames(output_basic))
expect_true("IsotopeLabelType" %in% colnames(output_basic))
expect_true("Condition"        %in% colnames(output_basic))
expect_true("BioReplicate"     %in% colnames(output_basic))

# Without heavyLabel all rows should be "L" (backwards compatible default)
expect_true(all(output_basic$IsotopeLabelType == "L"))

# Condition values should reflect R.Condition (0d, 8d, 32d)
expect_true(all(c("0d", "8d", "32d") %in% unique(output_basic$Condition)))

# BioReplicate falls back to R.Condition since R.Replicate is absent
expect_true(all(output_basic$BioReplicate %in% output_basic$Condition))

# Protein names come from PG.ProteinAccessions (no PG.ProteinGroups column)
expect_true(nrow(output_basic) > 0)
expect_false(any(is.na(output_basic$ProteinName)))


# --- Heavy label classification (Lys6) --------------------------------------

output_heavy = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "MS1Quantity",
    heavyLabel   = "Lys6",
    use_log_file = FALSE
)

expect_true("IsotopeLabelType" %in% colnames(output_heavy))

# Both heavy and light peptides must be present
expect_true("H" %in% unique(output_heavy$IsotopeLabelType))
expect_true("L" %in% unique(output_heavy$IsotopeLabelType))

# Heavy peptides must have [Lys6] in their PeptideSequence
heavy_rows = output_heavy[IsotopeLabelType == "H"]
expect_true(all(grepl("[Lys6]", heavy_rows$PeptideSequence, fixed = TRUE)))

# Light peptides must NOT have [Lys6] in their PeptideSequence
light_rows = output_heavy[IsotopeLabelType == "L"]
expect_false(any(grepl("[Lys6]", light_rows$PeptideSequence, fixed = TRUE)))

# Unlabeled (NA) peptides must NOT have [Lys6] in their PeptideSequence
na_rows = output_heavy[is.na(IsotopeLabelType)]
if (nrow(na_rows) > 0) {
    expect_false(any(grepl("[Lys6]", na_rows$PeptideSequence, fixed = TRUE)))
}


# --- MS2 intensity channel ---------------------------------------------------

output_ms2 = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "MS2Quantity",
    use_log_file = FALSE
)

expect_true("Intensity" %in% colnames(output_ms2))
# MS1 and MS2 intensity values should generally differ
output_ms1 = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "MS1Quantity",
    use_log_file = FALSE
)
# At least some intensities must differ between MS1 and MS2 channels
shared_keys = merge(
    output_ms1[, .(ProteinName, PeptideSequence, Run, ms1 = Intensity)],
    output_ms2[, .(ProteinName, PeptideSequence, Run, ms2 = Intensity)],
    by = c("ProteinName", "PeptideSequence", "Run")
)
if (nrow(shared_keys) > 0) {
    expect_false(all(shared_keys$ms1 == shared_keys$ms2, na.rm = TRUE))
}


# --- Raw column name as intensity string -------------------------------------

output_raw_col = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "FG.MS1Quantity",
    use_log_file = FALSE
)
# Should produce the same intensities as the alias "MS1Quantity"
expect_equal(nrow(output_raw_col), nrow(output_ms1))


# --- Invalid intensity column ------------------------------------------------

expect_error(
    SpectronauttoMSstatsFormat(
        boxcar_raw,
        intensity    = "FG.NonExistentColumn",
        use_log_file = FALSE
    ),
    "not found in input data"
)


# --- Novel heavy labels (non-Lys6) ------------------------------------------

# Deuterium leucine label; none in this file, so all should be L or NA
output_leu = SpectronauttoMSstatsFormat(
    boxcar_raw,
    intensity    = "MS1Quantity",
    heavyLabel   = "Leu6",
    use_log_file = FALSE
)
expect_false("H" %in% unique(output_leu$IsotopeLabelType))


# --- Backwards compatibility: standard Spectronaut file ----------------------
# The standard format (with all columns present) must still work unchanged.

spectronaut_std_path = system.file(
    "tinytest/raw_data/Spectronaut/spectronaut_input.csv",
    package = "MSstatsConvert")
spectronaut_std = data.table::fread(spectronaut_std_path)

output_std = SpectronauttoMSstatsFormat(spectronaut_std, use_log_file = FALSE)
expect_equal(ncol(output_std), 11)
expect_equal(nrow(output_std), 372)
expect_true(all(output_std$IsotopeLabelType == "L"))
