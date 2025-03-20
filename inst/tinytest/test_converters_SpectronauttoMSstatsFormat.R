# Test SpectronauttoMSstatsFormat ---------------------------
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
output = SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 372)
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


# Test SpectronauttoMSstatsFormat Missing Columns ---------------------------
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$F.ExcludedFromQuantification = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: FExcludedFromQuantification"
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$F.FrgLossType = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: FFrgLossType"
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$PG.ProteinGroups = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: PGProteinGroups"
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$FG.Charge = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: FGCharge"
)

