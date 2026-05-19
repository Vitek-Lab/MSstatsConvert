# Test MZMinetoMSstatsFormat ---------------------------
input_file_path = system.file("tinytest/raw_data/MZMine/mzmine_input.csv",
                              package = "MSstatsConvert")
annotation_file_path = system.file("tinytest/raw_data/MZMine/annotation.csv",
                                   package = "MSstatsConvert")
mzmine_ann_file_path = system.file("tinytest/raw_data/MZMine/mzmine_annotations.csv",
                                   package = "MSstatsConvert")
input = data.table::fread(input_file_path)
annot = data.table::fread(annotation_file_path)
mzmine_ann = data.table::fread(mzmine_ann_file_path)

# With mzmine_annotations supplied -------------------------------------------
output = MZMinetoMSstatsFormat(input, annotation = annot,
                               mzmine_annotations = mzmine_ann,
                               use_log_file = FALSE)
output_dt = data.table::as.data.table(output)

# Basic structure: 6 features x 4 runs = 24 rows, 11 standard columns
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 24)
expect_true("Run" %in% colnames(output))
expect_true("ProteinName" %in% colnames(output))
expect_true("PeptideSequence" %in% colnames(output))
expect_true("PrecursorCharge" %in% colnames(output))
expect_true("Intensity" %in% colnames(output))
expect_true("FragmentIon" %in% colnames(output))
expect_true("ProductCharge" %in% colnames(output))
expect_true("IsotopeLabelType" %in% colnames(output))
expect_true("Condition" %in% colnames(output))
expect_true("BioReplicate" %in% colnames(output))
expect_true("Fraction" %in% colnames(output))

# Metabolomics has no isotope labeling, so every row is "Light"
expect_true(all(output_dt$IsotopeLabelType == "Light"))

# Charge / fragment columns are not applicable for metabolomics
expect_true(all(is.na(output_dt$PrecursorCharge)))
expect_true(all(is.na(output_dt$FragmentIon)))
expect_true(all(is.na(output_dt$ProductCharge)))

# Fraction filled to 1
expect_true(all(output_dt$Fraction == 1))

# Annotation join: feature 2 has two annotation rows; the highest-scoring one wins
feature2_proteins = unique(output_dt[PeptideSequence == "2", ProteinName])
expect_equal(as.character(feature2_proteins), "GlucoseHigh")

# Clean annotation cases
feature1_proteins = unique(output_dt[PeptideSequence == "1", ProteinName])
expect_equal(as.character(feature1_proteins), "Caffeine")
feature3_proteins = unique(output_dt[PeptideSequence == "3", ProteinName])
expect_equal(as.character(feature3_proteins), "Lactate")
feature6_proteins = unique(output_dt[PeptideSequence == "6", ProteinName])
expect_equal(as.character(feature6_proteins), "Caffeine")

# Features without annotation rows fall back to the mz_rt string
feature4_proteins = unique(output_dt[PeptideSequence == "4", ProteinName])
expect_equal(as.character(feature4_proteins), "489.334_7.89")
feature5_proteins = unique(output_dt[PeptideSequence == "5", ProteinName])
expect_equal(as.character(feature5_proteins), "555.447_9.1")

# Zero-intensity input cells are converted to NA in output
# Feature 3 sampleB = 0  ->  NA
feature3_sampleB_int = output_dt[PeptideSequence == "3" & Run == "sampleBmzML",
                                  Intensity]
expect_true(is.na(feature3_sampleB_int))
# Feature 5 sampleB/C/D all = 0  ->  NA
feature5_zero_ints = output_dt[PeptideSequence == "5" &
                                   Run %in% c("sampleBmzML", "sampleCmzML", "sampleDmzML"),
                               Intensity]
expect_true(all(is.na(feature5_zero_ints)))

# Annotation merges correctly: sampleA is Control rep 1
sampleA_cond = unique(output_dt[Run == "sampleAmzML", Condition])
expect_equal(as.character(sampleA_cond), "Control")
sampleA_rep = unique(output_dt[Run == "sampleAmzML", BioReplicate])
expect_equal(as.character(sampleA_rep), "1")
sampleC_cond = unique(output_dt[Run == "sampleCmzML", Condition])
expect_equal(as.character(sampleC_cond), "Treatment")

# Intensity values trace back to input
feature1_sampleA_int = output_dt[PeptideSequence == "1" & Run == "sampleAmzML",
                                  Intensity]
expect_equal(as.numeric(feature1_sampleA_int), 1000)
feature2_sampleC_int = output_dt[PeptideSequence == "2" & Run == "sampleCmzML",
                                  Intensity]
expect_equal(as.numeric(feature2_sampleC_int), 5200)

# Without mzmine_annotations -------------------------------------------------
output_nolib = MZMinetoMSstatsFormat(input, annotation = annot,
                                     mzmine_annotations = NULL,
                                     use_log_file = FALSE)
output_nolib_dt = data.table::as.data.table(output_nolib)

# Every ProteinName is the mz_rt fallback string
expect_equal(ncol(output_nolib), 11)
expect_equal(nrow(output_nolib), 24)
expected_mz_rt = c("123.056_1.23", "245.129_3.45", "367.201_5.67",
                   "489.334_7.89", "555.447_9.1", "123.056_1.45")
expect_true(all(as.character(output_nolib_dt$ProteinName) %in% expected_mz_rt))
# Compound names from the library must not leak in
expect_false(any(as.character(output_nolib_dt$ProteinName) %in%
                 c("Caffeine", "GlucoseHigh", "GlucoseLow", "Lactate")))

# removeProtein_with1Feature filters non-Caffeine proteins -------------------
# Caffeine has 2 features (PeptideSequence "1" and "6"); all others have 1.
output_filtered = MZMinetoMSstatsFormat(input, annotation = annot,
                                        mzmine_annotations = mzmine_ann,
                                        removeProtein_with1Feature = TRUE,
                                        use_log_file = FALSE)
output_filtered_dt = data.table::as.data.table(output_filtered)

expect_equal(unique(as.character(output_filtered_dt$ProteinName)), "Caffeine")
# 2 features x 4 runs = 8 rows
expect_equal(nrow(output_filtered), 8)
expect_equal(sort(unique(as.character(output_filtered_dt$PeptideSequence))),
             c("1", "6"))
