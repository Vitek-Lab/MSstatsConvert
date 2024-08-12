.testHappyPath = function(input) {
    output = MSstatsConvert:::.cleanRawProteinProspector(input)
    expect_equal(ncol(output), 7)
    expect_equal(nrow(output), 966)
    expected_column_names = c(
        "ProteinName",
        "PeptideSequence",
        "PrecursorCharge",
        "Run",
        "PSM",
        "Channel",
        "Intensity"
    )
    missing_columns = setdiff(expected_column_names, colnames(output))
    expect_equal(length(missing_columns), 0)
}

# Test PD .cleanRawPDTMT ---------------------------

file_path = system.file("tinytest/processed_data/ProtProspector/MSstatsProtProspectorFilesObjectTMT.rds", package="MSstatsConvert")
input = readRDS(file_path)
.testHappyPath(input)

