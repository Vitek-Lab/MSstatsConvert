# Test SagetoMSstatsFormat ---------------------------
sage_path = system.file("tinytest/raw_data/Sage/lfq.tsv",
                        package = "MSstatsConvert")
annot_path = system.file("tinytest/raw_data/Sage/annotation.csv",
                         package = "MSstatsConvert")
if (!nzchar(sage_path) || !nzchar(annot_path)) {
    exit_file("Sage fixtures not present in inst/tinytest/raw_data/Sage/")
}

sage_raw = data.table::fread(sage_path)
annotation = data.table::fread(annot_path)

fixed_cols = c("peptide", "charge", "proteins", "q_value", "score",
               "spectral_angle")
intensity_cols = setdiff(colnames(sage_raw), fixed_cols)

output = SagetoMSstatsFormat(sage_raw, annotation, use_log_file = FALSE)
output = data.table::as.data.table(output)

# The 11 expected MSstats columns exist
expected_cols = c("ProteinName", "PeptideSequence", "PrecursorCharge",
                  "FragmentIon", "ProductCharge", "IsotopeLabelType",
                  "Condition", "BioReplicate", "Run", "Fraction", "Intensity")
for (col in expected_cols) {
    expect_true(col %in% colnames(output))
}
expect_equal(ncol(output), 11L)

# Wide-to-long melt produced one Run per intensity column (8 columns -> 8 runs)
expect_equal(length(intensity_cols), 8L)
expect_equal(data.table::uniqueN(output$Run), length(intensity_cols))

# Run values map to the correct Condition and BioReplicate.
# Output Run is the standardized form, so standardize the annotation Run too.
annotation_std = data.table::as.data.table(annotation)
annotation_std[, Run := MSstatsConvert:::.standardizeColnames(Run)]
run_map = merge(
    unique(output[, list(Run, Condition, BioReplicate)]),
    annotation_std[, list(Run, Condition, BioReplicate)],
    by = "Run", suffixes = c("_out", "_annot"))
expect_equal(nrow(run_map), length(intensity_cols))
expect_true(all(as.character(run_map$Condition_out) ==
                    as.character(run_map$Condition_annot)))
expect_true(all(as.character(run_map$BioReplicate_out) ==
                    as.character(run_map$BioReplicate_annot)))

# Zero intensities became NA, not 0
expect_false(any(output$Intensity == 0, na.rm = TRUE))
expect_true(any(is.na(output$Intensity)))

# IsotopeLabelType is "L" everywhere, and the annotation's own IsotopeLabelType
# column neither broke the merge nor produced a duplicate column
expect_true("IsotopeLabelType" %in% colnames(annotation))
expect_true(all(output$IsotopeLabelType == "L"))
expect_equal(sum(colnames(output) == "IsotopeLabelType"), 1L)

# FragmentIon and ProductCharge are NA
expect_true(all(is.na(output$FragmentIon)))
expect_true(all(is.na(output$ProductCharge)))

# q-value filter is load-bearing. Isolate its effect with removeFewMeasurements
# = FALSE so the only difference between the two runs is filter_with_Qvalue.
high_q_peptides = sage_raw[q_value > 0.01, unique(peptide)]
expect_true(length(high_q_peptides) >= 1L)

out_filter = SagetoMSstatsFormat(sage_raw, annotation,
                                 filter_with_Qvalue = TRUE,
                                 removeFewMeasurements = FALSE,
                                 use_log_file = FALSE)
out_filter = data.table::as.data.table(out_filter)
out_nofilter = SagetoMSstatsFormat(sage_raw, annotation,
                                   filter_with_Qvalue = FALSE,
                                   removeFewMeasurements = FALSE,
                                   use_log_file = FALSE)
out_nofilter = data.table::as.data.table(out_nofilter)

signal_with_filter = unique(out_filter[!is.na(Intensity), PeptideSequence])
signal_without_filter = unique(out_nofilter[!is.na(Intensity), PeptideSequence])

# With filter FALSE the 5 high-q peptides keep signal; with TRUE they are gone
expect_true(all(high_q_peptides %in% signal_without_filter))
expect_false(any(high_q_peptides %in% signal_with_filter))

# The peptides that lose all signal when the filter is switched on are exactly
# the high-q peptides
expect_equal(sort(setdiff(signal_without_filter, signal_with_filter)),
             sort(high_q_peptides))


# Test SagetoMSstatsFormat on the charge-resolved fixture -------------------
sage_cr_path = system.file("tinytest/raw_data/Sage/lfq_charge_resolved.tsv",
                           package = "MSstatsConvert")
annot_cr_path = system.file("tinytest/raw_data/Sage/annotation_charge_resolved.csv",
                            package = "MSstatsConvert")
if (nzchar(sage_cr_path) && nzchar(annot_cr_path)) {
    sage_cr = data.table::fread(sage_cr_path)
    annot_cr = data.table::fread(annot_cr_path)

    # Single run: keep single-measurement features so nothing is dropped
    out_cr = SagetoMSstatsFormat(sage_cr, annot_cr,
                                 removeFewMeasurements = FALSE,
                                 use_log_file = FALSE)
    out_cr = data.table::as.data.table(out_cr)

    # Real charges come through, and -1 does not appear
    expect_true(all(c(2, 3) %in% unique(out_cr$PrecursorCharge)))
    expect_false(-1 %in% unique(out_cr$PrecursorCharge))

    # A peptide seen at two charges yields two distinct features
    charge_counts = out_cr[, list(n_charge = data.table::uniqueN(PrecursorCharge)),
                           by = PeptideSequence]
    expect_true(any(charge_counts$n_charge == 2L))
    two_charge_pep = charge_counts[n_charge == 2L, PeptideSequence][1]
    feats = unique(out_cr[PeptideSequence == two_charge_pep,
                          list(PeptideSequence, PrecursorCharge)])
    expect_equal(nrow(feats), 2L)
}
