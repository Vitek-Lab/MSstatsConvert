# Test DIANN .cleanRawDIANN ---------------------------

.validateOutput = function(output) {
    expect_equal(ncol(output), 12)
    expect_true(nrow(output) > 0)
    expect_true("Run" %in% colnames(output))
    expect_true("ProteinName" %in% colnames(output))
    expect_true("PeptideSequence" %in% colnames(output))
    expect_true("PrecursorCharge" %in% colnames(output))
    expect_true("Intensity" %in% colnames(output))
    expect_true("DetectionQValue" %in% colnames(output))
    expect_true("PrecursorMz" %in% colnames(output))
    expect_true("LibQValue" %in% colnames(output))
    expect_true("LibPGQValue" %in% colnames(output))
    expect_true("FragmentInfo" %in% colnames(output))
    expect_true("FragmentIon" %in% colnames(output))
    expect_true("ProductCharge" %in% colnames(output))
    expect_equal(class(output$Intensity), "numeric")
}

file_path = system.file("tinytest/processed_data/DIANN/MSstatsDIANNFilesObject.rds", package="MSstatsConvert")
input = readRDS(file_path)
output = MSstatsConvert:::.cleanRawDIANN(input)
.validateOutput(output)
output = MSstatsConvert:::.cleanRawDIANN(input, quantificationColumn = "FragmentQuantRaw")
.validateOutput(output)

