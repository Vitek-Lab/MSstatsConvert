# Test DIANNtoMSstatsFormat SILAC protein turnover ---------------------------
# Dataset source: MSV000097050
input_file_path_silac = system.file("tinytest/raw_data/DIANN/diann_input_silac.csv",
                                    package = "MSstatsConvert")
annotation_file_path_silac = system.file("tinytest/raw_data/DIANN/annotation_silac.csv",
                                         package = "MSstatsConvert")
input_silac = data.table::fread(input_file_path_silac)
annot_silac = data.table::fread(annotation_file_path_silac)
output_silac = DIANNtoMSstatsFormat(input_silac, annotation = annot_silac,
                                    labeledAminoAcids = c("K"),
                                    removeFewMeasurements = FALSE,
                                    use_log_file = FALSE)
output_silac_dt = data.table::as.data.table(output_silac)

# Basic structure
expect_equal(ncol(output_silac), 11)
expect_equal(nrow(output_silac), 1020)
expect_true("Run" %in% colnames(output_silac))
expect_true("ProteinName" %in% colnames(output_silac))
expect_true("PeptideSequence" %in% colnames(output_silac))
expect_true("PrecursorCharge" %in% colnames(output_silac))
expect_true("Intensity" %in% colnames(output_silac))
expect_true("FragmentIon" %in% colnames(output_silac))
expect_true("ProductCharge" %in% colnames(output_silac))
expect_true("IsotopeLabelType" %in% colnames(output_silac))
expect_true("Condition" %in% colnames(output_silac))
expect_true("BioReplicate" %in% colnames(output_silac))
expect_true("Fraction" %in% colnames(output_silac))

# IsotopeLabelType is classified as H/L/NA — not the "Light" default from
# the labeledAminoAcids=NULL path
expect_false("Light" %in% output_silac$IsotopeLabelType)
label_counts = table(output_silac$IsotopeLabelType, useNA = "ifany")
expect_equal(unname(label_counts["H"]), 350L)
expect_equal(unname(label_counts["L"]), 350L)
expect_equal(sum(is.na(output_silac$IsotopeLabelType)), 320L)

# SILAC modification tags are stripped from PeptideSequence
expect_false(any(grepl("SILAC", output_silac$PeptideSequence)))

# Peptides without a labeled K get NA IsotopeLabelType
unlabeled_rows = output_silac_dt[PeptideSequence %in% c("AAAAADLANR", "AAAADGEPLHNEEER")]
expect_true(all(is.na(unlabeled_rows$IsotopeLabelType)))

# Peptides with a labeled K get only H or L
k_pep_rows = output_silac_dt[PeptideSequence %in% c("AAAAAAAAQMHTK", "AAAAAAAK",
                                                     "AAAAAAK", "AAAAAGGK")]
expect_true(all(k_pep_rows$IsotopeLabelType %in% c("H", "L")))

# Annotation merges correctly: run dAL_AT_Long_B1 maps to Condition 10
cond_b1 = unique(output_silac_dt[Run == "dAL_AT_Long_B1", Condition])
expect_equal(as.character(cond_b1), "10")

# Heavy fragment intensities trace back to the heavy input rows.
# Input: AAAAAAAAQMHTK(SILAC-K-H), Run dAL_AT_Long_B1 has
#   Fragment.Quant.Corrected = "0;18378.38477;17024.44141;15426.27441;38014.78906;"
#   Fragment.Info ions (after NH3 removal):
#     y8^1/863.4499512 -> 18378.38477, y9^1/934.4870605 -> 17024.44141,
#     y6^1/721.3757324 -> 15426.27441
heavy_b1 = output_silac_dt[PeptideSequence == "AAAAAAAAQMHTK" &
                            IsotopeLabelType == "H" &
                            Run == "dAL_AT_Long_B1"]
heavy_y8 = heavy_b1[FragmentIon == "y8^1/863.4499512", Intensity]
expect_equal(heavy_y8, 18378.38477, tolerance = 1)
heavy_y9 = heavy_b1[FragmentIon == "y9^1/934.4870605", Intensity]
expect_equal(heavy_y9, 17024.44141, tolerance = 1)
heavy_y6 = heavy_b1[FragmentIon == "y6^1/721.3757324", Intensity]
expect_equal(heavy_y6, 15426.27441, tolerance = 1)

# Light fragment intensities trace back to the light input rows.
# Input: AAAAAAAAQMHTK(SILAC-K-L), Run dAL_AT_Long_B1 has
#   Fragment.Quant.Corrected = "120170.8125;39879.93359;13874.50293;33449.03516;0;33788.10156;"
#   Fragment.Info ions (after NH3 removal):
#     y5^1/644.3184814 -> 120170.8125, y8^1/857.4298096 -> 39879.93359,
#     y10^1/999.5040283 -> 13874.50293, y9^1/928.4669189 -> 33449.03516
light_b1 = output_silac_dt[PeptideSequence == "AAAAAAAAQMHTK" &
                            IsotopeLabelType == "L" &
                            Run == "dAL_AT_Long_B1"]
light_y5 = light_b1[FragmentIon == "y5^1/644.3184814", Intensity]
expect_equal(light_y5, 120170.8125, tolerance = 1)
light_y8 = light_b1[FragmentIon == "y8^1/857.4298096", Intensity]
expect_equal(light_y8, 39879.93359, tolerance = 1)
light_y10 = light_b1[FragmentIon == "y10^1/999.5040283", Intensity]
expect_equal(light_y10, 13874.50293, tolerance = 1)
light_y9 = light_b1[FragmentIon == "y9^1/928.4669189", Intensity]
expect_equal(light_y9, 33449.03516, tolerance = 1)

# Zero-intensity input fragments are converted to NA in output
# AAAAAAAAQMHTK(SILAC-K-L) B1 has y6^1/715.3555908 = 0 in input
light_y6_zero = light_b1[FragmentIon == "y6^1/715.3555908", Intensity]
expect_true(is.na(light_y6_zero))

# Verify cross-run consistency: B8 light intensities also trace to input
# Input: AAAAAAAAQMHTK(SILAC-K-L), Run dAL_AT_Long_B8 has
#   Fragment.Quant.Corrected = "27676.55078;5134.108887;4891.382324;3844.573975;3443.866699;8253.693359;"
light_b8 = output_silac_dt[PeptideSequence == "AAAAAAAAQMHTK" &
                            IsotopeLabelType == "L" &
                            Run == "dAL_AT_Long_B8"]
light_y5_b8 = light_b8[FragmentIon == "y5^1/644.3184814", Intensity]
expect_equal(light_y5_b8, 27676.55078, tolerance = 1)
light_y8_b8 = light_b8[FragmentIon == "y8^1/857.4298096", Intensity]
expect_equal(light_y8_b8, 5134.108887, tolerance = 1)

# Test DIANNtoMSstatsFormat PXD055942 (DIANN 2.0 + Channel-based SILAC) ------
# Dataset: PXD055942 (heart protein turnover, 0/8/32 day timepoints).
# Uses quantificationColumn = "auto" (Fr.N.Quantity columns) and the Channel
# column path for IsotopeLabelType (no SILAC suffix in ModifiedSequence).
input_file_path_pxd = system.file("tinytest/raw_data/DIANN/PXD055942.parquet",
                                   package = "MSstatsConvert")
annotation_file_path_pxd = system.file("tinytest/raw_data/DIANN/annotation_PXD055942.csv",
                                        package = "MSstatsConvert")
input_pxd = arrow::read_parquet(input_file_path_pxd)
annot_pxd = data.table::fread(annotation_file_path_pxd)
output_pxd = DIANNtoMSstatsFormat(input_pxd, annotation = annot_pxd,
                                   labeledAminoAcids = c("K"),
                                   removeFewMeasurements = FALSE,
                                   quantificationColumn = "auto",
                                   use_log_file = FALSE)
output_pxd_dt = data.table::as.data.table(output_pxd)

# Basic structure
expect_equal(ncol(output_pxd), 11)
expect_equal(nrow(output_pxd), 20400)
expect_true("Run" %in% colnames(output_pxd))
expect_true("ProteinName" %in% colnames(output_pxd))
expect_true("PeptideSequence" %in% colnames(output_pxd))
expect_true("PrecursorCharge" %in% colnames(output_pxd))
expect_true("Intensity" %in% colnames(output_pxd))
expect_true("FragmentIon" %in% colnames(output_pxd))
expect_true("ProductCharge" %in% colnames(output_pxd))
expect_true("IsotopeLabelType" %in% colnames(output_pxd))
expect_true("Condition" %in% colnames(output_pxd))
expect_true("BioReplicate" %in% colnames(output_pxd))
expect_true("Fraction" %in% colnames(output_pxd))

# Channel path used: H/L/NA counts, no "Light" default
expect_false("Light" %in% output_pxd$IsotopeLabelType)
pxd_label_counts = table(output_pxd$IsotopeLabelType, useNA = "ifany")
expect_equal(unname(pxd_label_counts["H"]), 8712L)
expect_equal(unname(pxd_label_counts["L"]), 8712L)
expect_equal(sum(is.na(output_pxd$IsotopeLabelType)), 2976L)

# Channel path does NOT strip (SILAC) from PeptideSequence — that notation is
# part of the peptide identity in this format, not a label indicator
expect_true(any(grepl("(SILAC)", output_pxd$PeptideSequence, fixed = TRUE)))

# Unlabeled peptides (no K, no SILAC) get NA IsotopeLabelType
pxd_unlabeled = output_pxd_dt[PeptideSequence %in% c("AVLEEAEFQR", "DDEGLYTLR",
                                                       "DTELAEELLQWFLQEEK")]
expect_true(all(is.na(pxd_unlabeled$IsotopeLabelType)))

# Annotation maps correctly to all three timepoint conditions
expect_true("0 day" %in% output_pxd$Condition)
expect_true("8 days" %in% output_pxd$Condition)
expect_true("32 days" %in% output_pxd$Condition)
run_0d = "20210805_BoxCarmax1st_wideMS1_JM_pSIL_pro_heart_0d_1"
expect_equal(as.character(unique(output_pxd_dt[Run == run_0d, Condition])), "0 day")

# Heavy fragment intensities match Fr.N.Quantity from the Channel = "H" input row.
# Input: AFMTADLPNELIELLEK(SILAC)(SILAC)(SILAC), Run 32d_1, Channel H has
#   Fr.0.Quantity = 662450.1 (largest fragment — should be the max output intensity)
h_pxd_pep = "AFMTADLPNELIELLEK(SILAC)(SILAC)(SILAC)"
h_pxd_run = "20210805_BoxCarmax2nd_wideMS1_JM_pSIL_pro_heart_32d_1"
h_pxd_ints = output_pxd_dt[PeptideSequence == h_pxd_pep &
                             IsotopeLabelType == "H" &
                             Run == h_pxd_run, Intensity]
expect_equal(max(h_pxd_ints, na.rm = TRUE), 662450.1, tolerance = 1)

# Light fragment intensities match Fr.N.Quantity from the Channel = "L" input row.
# Input: APVYSGSSPVSGYFVDFK(SILAC)(SILAC)(SILAC)EEDSGEWK(SILAC)(SILAC)(SILAC),
#   Run 32d_2, Channel L has Fr.1.Quantity = 116961.9 (largest fragment)
l_pxd_pep = "APVYSGSSPVSGYFVDFK(SILAC)(SILAC)(SILAC)EEDSGEWK(SILAC)(SILAC)(SILAC)"
l_pxd_run = "20210805_BoxCarmax1st_wideMS1_JM_pSIL_pro_heart_32d_2"
l_pxd_ints = output_pxd_dt[PeptideSequence == l_pxd_pep &
                             IsotopeLabelType == "L" &
                             Run == l_pxd_run, Intensity]
expect_equal(max(l_pxd_ints, na.rm = TRUE), 116961.9, tolerance = 1)
# Fr.0.Quantity = 4973.411 also appears in the output
expect_true(any(abs(l_pxd_ints - 4973.411) < 1, na.rm = TRUE))

# Test DIANNtoMSstatsFormat ---------------------------
input_file_path = system.file("tinytest/raw_data/DIANN/diann_input.tsv", package="MSstatsConvert")
annotation_file_path = system.file("tinytest/raw_data/DIANN/annotation.csv", package = "MSstatsConvert")
input = data.table::fread(input_file_path)
annot = data.table::fread(annotation_file_path)
output = DIANNtoMSstatsFormat(input, annotation = annot, MBR = FALSE, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 348) 
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
# When labeledAminoAcids = NULL (default), IsotopeLabelType is filled with "Light"
expect_true(all(output$IsotopeLabelType == "Light"))

# Test DIANNtoMSstatsFormat DIANN 2.0 ------------------------
input_file_path = system.file("tinytest/raw_data/DIANN/diann_2.0.parquet", package="MSstatsConvert")
annotation_file_path = system.file("tinytest/raw_data/DIANN/annotation_diann_2.0.csv", package = "MSstatsConvert")
input = arrow::read_parquet(input_file_path)
annot = data.table::fread(annotation_file_path)
output = DIANNtoMSstatsFormat(input, annotation = annot, MBR = FALSE, use_log_file = FALSE, quantificationColumn = 'auto')
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 192)
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
# When labeledAminoAcids = NULL (default), IsotopeLabelType is filled with "Light"
expect_true(all(output$IsotopeLabelType == "Light"))