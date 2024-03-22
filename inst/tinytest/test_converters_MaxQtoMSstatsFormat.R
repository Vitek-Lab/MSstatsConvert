# Test MaxQtoMSstatsFormat ---------------------------
mq_ev = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_ev.csv",
                                      package = "MSstatsConvert"))
mq_pg = data.table::fread(system.file("tinytest/raw_data/MaxQuant/mq_pg.csv",
                                      package = "MSstatsConvert"))
annot = data.table::fread(system.file("tinytest/raw_data/MaxQuant/annotation.csv",
                                      package = "MSstatsConvert"))
output = MaxQtoMSstatsFormat(mq_ev, annot, mq_pg, use_log_file = FALSE)
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