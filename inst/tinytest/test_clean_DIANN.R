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

# Q-value filtering
expect_qvalue_cutoff <- function(output, col, cutoff) {
    expect_equal(
        sum(output[[col]] > cutoff),
        sum(output[["Intensity"]] == 0 & output[[col]] > cutoff),
        info = sprintf(
            "All rows with %s > %s should have %s == 0",
            col, cutoff, "Intensity"
        )
    )
    expect_equal(
        sum(output[[col]] <= cutoff),
        nrow(output) - sum(output[[col]] > cutoff),
        info = sprintf(
            "Rows with %s <= %s should account for all rows not above the cutoff",
            col, cutoff
        )
    )
}
output <- MSstatsConvert:::.cleanRawDIANN(input, global_qvalue_cutoff = 0.005)
expect_qvalue_cutoff(output, "DetectionQValue", 0.005)
output <- MSstatsConvert:::.cleanRawDIANN(input, qvalue_cutoff = 0.00001)
expect_qvalue_cutoff(output, "LibQValue", 0.00001)
output <- MSstatsConvert:::.cleanRawDIANN(input, pg_qvalue_cutoff = 0.001)
expect_qvalue_cutoff(output, "LibPGQValue", 0.001)
output <- MSstatsConvert:::.cleanRawDIANN(input, MBR = FALSE, qvalue_cutoff = 0.001)
expect_qvalue_cutoff(output, "GlobalQValue", 0.001)
output <- MSstatsConvert:::.cleanRawDIANN(input, MBR = FALSE, pg_qvalue_cutoff = 0.0002)
expect_qvalue_cutoff(output, "GlobalPGQValue", 0.0002)

# Test .assignDIANNIsotopeLabelType ---------------------------

dt_channel = data.table::data.table(
    PeptideSequence = c("PEPTIDEK", "PEPTIDEK", "PEPTIDEK"),
    Channel = c("H", "L", "other")
)
result_channel = MSstatsConvert:::.assignDIANNIsotopeLabelType(
    dt_channel, labeledAminoAcids = c("K"), has_channel = TRUE
)
expect_equal(result_channel$IsotopeLabelType, c("H", "L", NA_character_))
expect_false("Channel" %in% colnames(result_channel))

dt_null = data.table::data.table(
    PeptideSequence = c("PEPTIDEK(SILAC-K-H)", "PEPTIDEK(SILAC-K-L)")
)
result_null = MSstatsConvert:::.assignDIANNIsotopeLabelType(
    dt_null, labeledAminoAcids = NULL, has_channel = FALSE
)
expect_false("IsotopeLabelType" %in% colnames(result_null))
expect_equal(nrow(result_null), 2L)

dt_label = data.table::data.table(
    PeptideSequence = c("PEPTIDEK(SILAC-K-H)", "PEPTIDEK(SILAC-K-L)",
                        "PEPTIDEK(label-K-H)", "PEPTIDEK(label-K-L)",
                        "PEPTIDEK")
)
result_label = MSstatsConvert:::.assignDIANNIsotopeLabelType(
    dt_label, labeledAminoAcids = c("K"), has_channel = FALSE
)
expect_equal(result_label$IsotopeLabelType, c("H", "L", "H", "L", NA_character_))
expect_equal(unique(result_label$PeptideSequence), "PEPTIDEK")

dt_multi_aa = data.table::data.table(
    PeptideSequence = c("PEPTIDEK(SILAC-K-H)", "PEPTIDER(SILAC-R-H)",
                        "PEPTIDEK(SILAC-K-L)", "PEPTIDER(SILAC-R-L)",
                        "PEPTIDEAC")
)
result_multi_aa = MSstatsConvert:::.assignDIANNIsotopeLabelType(
    dt_multi_aa, labeledAminoAcids = c("K", "R"), has_channel = FALSE
)
expect_equal(result_multi_aa$IsotopeLabelType,
             c("H", "H", "L", "L", NA_character_))
expect_equal(sort(unique(result_multi_aa$PeptideSequence)),
             c("PEPTIDEAC", "PEPTIDEK", "PEPTIDER"))
