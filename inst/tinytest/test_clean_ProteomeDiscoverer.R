.testHappyPath = function(input) {
    output = MSstatsConvert:::.cleanRawPDTMT(input)
    expect_equal(ncol(output), 11)
    expect_true(nrow(output) > 0)
    expected_column_names = c(
        "ProteinName",
        "numProtein",
        "PeptideSequence",
        "PrecursorCharge",
        "IonsScore",
        "Run",
        "QuanInfo",
        "IsolationInterference",
        "PSM",
        "Channel",
        "Intensity"
    )
    missing_columns = setdiff(expected_column_names, colnames(output))
    expect_equal(length(missing_columns), 0)
}

# Test PD .cleanRawPDTMT ---------------------------

file_path = system.file("tinytest/processed_data/PD/MSstatsProteomeDiscovererFilesObjectTMT.rds", package="MSstatsConvert")
input = readRDS(file_path)
.testHappyPath(input)

# Test PD .cleanRawPDTMT Column Validation Error ---------------------------
file_path = system.file("tinytest/processed_data/PD/MSstatsProteomeDiscovererFilesObjectTMTBadInput.rds", 
                        package="MSstatsConvert")
input_with_missing_column = readRDS(file_path)
expect_error(MSstatsConvert:::.cleanRawPDTMT(input_with_missing_column),
             pattern = "The following columns are missing from the input data")

# Test PD .cleanRawPDTMT for PD 3.1+, where "# Proteins" & "# Protein Groups" 
# columns are renamed to "Number of Proteins" & "Number of Protein Groups"
file_path = system.file("tinytest/processed_data/PD/MSstatsProteomeDiscovererFilesObjectTMT3-1Format.rds", package="MSstatsConvert")
input = readRDS(file_path)
.testHappyPath(input)

