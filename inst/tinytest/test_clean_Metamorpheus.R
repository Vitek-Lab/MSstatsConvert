# Test Metamorpheus .cleanRawMetamorpheus ---------------------------

file_path = system.file("tinytest/processed_data/Metamorpheus/MSstatsMetamorpheusFilesObject.rds", package="MSstatsConvert")
input = readRDS(file_path)
output = MSstatsConvert:::.cleanRawMetamorpheus(input)
expect_equal(ncol(output), 5)
expect_true(nrow(output) > 0)
expect_true("Run" %in% colnames(output))
expect_true("ProteinName" %in% colnames(output))
expect_true("PeptideSequence" %in% colnames(output))
expect_true("PrecursorCharge" %in% colnames(output))
expect_true("Intensity" %in% colnames(output))