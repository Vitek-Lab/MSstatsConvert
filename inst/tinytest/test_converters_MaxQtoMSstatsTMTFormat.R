# Test MaxQtoMSstatsTMTFormat ---------------------------
mq_ev = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_ev.csv",
                                      package = "MSstatsConvert"))
mq_pg = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_pg.csv",
                                      package = "MSstatsConvert"))
annot = data.table::fread(system.file("tinytest/raw_data/MaxQuantTMT/mq_annotation.csv",
                                      package = "MSstatsConvert"))
output = MaxQtoMSstatsTMTFormat(mq_ev, mq_pg, annot, use_log_file = FALSE)
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

bad_annotation = data.table::fread(system.file("tinytest/raw_data/PDTMT/pd_annotation.csv",
                                              package = "MSstatsConvert"))

expect_error(MaxQtoMSstatsTMTFormat(evidence = mq_ev, proteinGroups = mq_pg,
                                    annotation = bad_annotation)) # wrong annotation file

# Verify output intensities are present in input
intensity_cols = grep("Reporter.intensity.corrected", colnames(mq_ev), value = TRUE)
input_intensities = as.character(unlist(mq_ev[, intensity_cols, with = FALSE]))
output_intensities = as.character(output$Intensity[!is.na(output$Intensity)])
expect_true(all(output_intensities %in% input_intensities))
