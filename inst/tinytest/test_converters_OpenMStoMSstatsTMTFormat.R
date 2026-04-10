# Test OpenMStoMSstatsTMTFormat ---------------------------
raw = data.table::fread(system.file("tinytest/raw_data/OpenMSTMT/openmstmt_input.csv",
                                    package = "MSstatsConvert"))
output = OpenMStoMSstatsTMTFormat(raw, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_true(nrow(output) > 0)
expect_true("Run" %in% colnames(output))
expect_true("ProteinName" %in% colnames(output))
expect_true("PeptideSequence" %in% colnames(output))
expect_true("Charge" %in% colnames(output))
expect_true("Intensity" %in% colnames(output))
expect_true("TechRepMixture" %in% colnames(output))
expect_true("PSM" %in% colnames(output))
expect_true("Mixture" %in% colnames(output))
expect_true("Condition" %in% colnames(output))
expect_true("BioReplicate" %in% colnames(output))
expect_true("Channel" %in% colnames(output))

# Verify output intensities are present in input
input_intensities = as.character(raw$Intensity)
output_intensities = as.character(output$Intensity[!is.na(output$Intensity)])
expect_true(all(output_intensities %in% input_intensities))
