# Test OpenSWATHtoMSstatsFormat ---------------------------
os_raw = system.file("tinytest/raw_data/OpenSWATH/openswath_input.csv",
                             package = "MSstatsConvert")
annot = system.file("tinytest/raw_data/OpenSWATH/annot_os.csv",
                    package = "MSstatsConvert")
os_raw = data.table::fread(os_raw)
annot = data.table::fread(annot)

output = OpenSWATHtoMSstatsFormat(os_raw, annot, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_true(nrow(output) > 0)
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