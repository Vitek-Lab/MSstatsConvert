# Test PhilosophertoMSstatsTMTFormat ---------------------------
input = data.table::fread(system.file("tinytest/raw_data/Philosopher/msstats.csv",
                                      package = "MSstatsConvert"))
annot = data.table::fread(system.file("tinytest/raw_data/Philosopher/MSstatsTMT_annotation.csv",
                                      package = "MSstatsConvert"))
output = PhilosophertoMSstatsTMTFormat(input, annot, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 550)
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
intensity_cols = grep("^Channel\\.", colnames(input), value = TRUE)
input_intensities = as.character(unlist(input[, intensity_cols, with = FALSE]))
output_intensities = as.character(output$Intensity[!is.na(output$Intensity)])
expect_true(all(output_intensities %in% input_intensities))
