# Test PDtoMSstatsTMTFormat ---------------------------
pd_raw = data.table::fread(system.file("tinytest/raw_data/PDTMT/pdtmt_input.csv",
                                       package = "MSstatsConvert"))
annot = data.table::fread(system.file("tinytest/raw_data/PDTMT/pd_annotation.csv",
                                      package = "MSstatsConvert"))
output = PDtoMSstatsTMTFormat(pd_raw, annot, use_log_file = FALSE)
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

raw.pd = data.table::fread(system.file("tinytest/raw_data/PDTMT/pdtmt_input.csv",
                                       package = "MSstatsConvert"))
annotation.pd = data.table::fread(system.file("tinytest/raw_data/PDTMT/pd_annotation.csv",
                                              package = "MSstatsConvert"))

expect_error(PDtoMSstatsTMTFormat(input = pd_raw[, !colnames(pd_raw) == "Protein.Accessions"], # missing columns in input
                                  annotation = annot, use_log_file = FALSE))

expect_error(PDtoMSstatsTMTFormat(input = pd_raw,
                                  annotation = annot[, !colnames(annot) == "Condition"], # missing columns in annotation
                                  use_log_file = FALSE))

# Verify output intensities are present in input
intensity_cols = grep("^Abundance", colnames(pd_raw), value = TRUE)
input_intensities = as.character(unlist(pd_raw[, intensity_cols, with = FALSE]))
output_intensities = as.character(output$Intensity[!is.na(output$Intensity)])
expect_true(all(output_intensities %in% input_intensities))
