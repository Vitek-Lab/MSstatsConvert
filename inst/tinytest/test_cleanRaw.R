# DIAUmpire
dia_frag = data.table::fread("./raw_data/DIAUmpire/dia_frag.csv")
dia_pept = data.table::fread("./raw_data/DIAUmpire/dia_pept.csv")
dia_prot = data.table::fread("./raw_data/DIAUmpire/dia_prot.csv")
dia_import = MSstatsConvert::MSstatsImport(list(Fragments = dia_frag,
                                                Peptides = dia_pept,
                                                Proteins = dia_prot),
                                           "MSstats", "DIAUmpire")
dia_import2 = MSstatsConvert::MSstatsImport(
    list(Fragments = dia_frag,
         Peptides = dia_pept[, colnames(dia_pept) != "Selected_fragments", 
                             with = FALSE],
         Proteins = dia_prot),
    "MSstats", "DIAUmpire")
dia_import3 = MSstatsConvert::MSstatsImport(
    list(Fragments = dia_frag,
         Peptides = dia_pept,
         Proteins = dia_prot[, colnames(dia_prot) != "Selected_peptides", 
                             with = FALSE]),
    "MSstats", "DIAUmpire")
dia_frag2 = data.table::copy(dia_frag)
dia_frag2$Fragment = gsub("\\+[0-9]+", "", dia_frag2$Fragment)
dia_import4 = MSstatsConvert::MSstatsImport(
    list(Fragments = dia_frag2,
         Peptides = dia_pept,
         Proteins = dia_prot),
    "MSstats", "DIAUmpire")
dia_cleaned = MSstatsConvert::MSstatsClean(dia_import, TRUE, TRUE)
dia_cleaned2 = MSstatsConvert::MSstatsClean(dia_import4, TRUE, TRUE)
dia_cleaned3 = MSstatsConvert::MSstatsClean(dia_import, TRUE, FALSE)
expect_error(MSstatsConvert::MSstatsClean(dia_import, FALSE, TRUE))
expect_equal(ncol(dia_cleaned), 5)
expect_equal(ncol(dia_cleaned3), 5)
expect_true(nrow(dia_cleaned) > 0)
expect_true(nrow(dia_cleaned3) > 0)
expect_error(MSstatsConvert::MSstatsClean(dia_import, FALSE, FALSE))
expect_error(MSstatsConvert::MSstatsClean(dia_import2, TRUE, TRUE))
expect_error(MSstatsConvert::MSstatsClean(dia_import3, TRUE, TRUE))
expect_equal(dia_cleaned, dia_cleaned2)
# MaxQuant
mq_ev = data.table::fread("./raw_data/MaxQuant/mq_ev.csv")
mq_pg = data.table::fread("./raw_data/MaxQuant/mq_pg.csv")
mq_import = MSstatsConvert::MSstatsImport(list(evidence = mq_ev, 
                                               protein_groups = mq_pg),
                                          "MSstats", "MaxQuant")
mq_cleaned = MSstatsConvert::MSstatsClean(mq_import, protein_id_col = "Proteins")
mq_cleaned_site = MSstatsConvert::MSstatsClean(mq_import, 
                                               protein_id_col = "Proteins",
                                               remove_by_site = TRUE)
expect_true(nrow(mq_cleaned_site) < nrow(mq_cleaned))
expect_true(nrow(mq_cleaned_site) > 0)
expect_equal(
    ncol(MSstatsConvert::MSstatsClean(mq_import, protein_id_col = "Proteins")),
    7
)
expect_true(nrow(mq_cleaned) > 0)
# MaxQuantTMT
mqtmt_ev = data.table::fread("./raw_data/MaxQuantTMT/mq_ev.csv")
mqtmt_pg = data.table::fread("./raw_data/MaxQuantTMT/mq_pg.csv")
mqtmt_import = MSstatsConvert::MSstatsImport(list(evidence = mqtmt_ev, 
                                                  protein_groups = mqtmt_pg),
                                             "MSstatsTMT", "MaxQuant")
mqtmt_cleaned = MSstatsConvert::MSstatsClean(mqtmt_import, 
                                             protein_id_col = "Proteins")
expect_equal(
    ncol(mqtmt_cleaned),
    8
)
expect_true(nrow(mq_cleaned) > 0)
# OpenMS
openms_input = data.table::fread("./raw_data/OpenMS/openms_input.csv")
openms_import = MSstatsConvert::MSstatsImport(list(input = openms_input), 
                                              "MSstats", "OpenMS")
openms_import2 = MSstatsConvert::MSstatsImport(list(input = openms_input[, -6]), 
                                               "MSstats", "OpenMS")
om_cleaned = MSstatsConvert::MSstatsClean(openms_import)
expect_equal(ncol(om_cleaned), 10)
expect_equal(ncol(MSstatsConvert::MSstatsClean(openms_import2)), 10)
expect_true(nrow(om_cleaned) > 0)
# OpenMSTMT
openmstmt_input = data.table::fread("./raw_data/OpenMSTMT/openmstmt_input.csv")
openmstmt_import = MSstatsConvert::MSstatsImport(list(input = openmstmt_input), 
                                                 "MSstatsTMT", "OpenMS")
omtmt_cleaned = MSstatsConvert::MSstatsClean(openmstmt_import)
expect_equal(ncol(omtmt_cleaned), 12)
expect_true(nrow(omtmt_cleaned) > 0)
# OpenSWATH
openswath_input = data.table::fread("./raw_data/OpenSWATH/openswath_input.csv")
openswath_import = MSstatsConvert::MSstatsImport(list(input = openswath_input), 
                                                 "MSstats", "OpenSWATH")
os_cleaned = MSstatsConvert::MSstatsClean(openswath_import)
expect_equal(ncol(os_cleaned), 8)
expect_true(nrow(os_cleaned) > 0)
# PD
pd_input = data.table::fread("./raw_data/PD/pd_input.csv")
pd_input_frac = data.table::copy(pd_input)
pd_input_frac$Fraction = 1
pd_import = MSstatsConvert::MSstatsImport(list(input = pd_input), 
                                          "MSstats", "ProteomeDiscoverer")
pd_import_frac = MSstatsConvert::MSstatsImport(list(input = pd_input_frac), 
                                               "MSstats", "ProteomeDiscoverer")
pd_cleaned = MSstatsConvert::MSstatsClean(
    pd_import, protein_id_column = "ProteinGroupAccessions",
    sequence_column = "Sequence", quantification_column = "Intensity",
    remove_shared = TRUE)
pd_cleaned_frac = MSstatsConvert::MSstatsClean(
    pd_import_frac, protein_id_column = "ProteinGroupAccessions",
    sequence_column = "Sequence", quantification_column = "Intensity",
    remove_shared = TRUE)
expect_equal(ncol(pd_cleaned), 5)
expect_true(nrow(pd_cleaned) > 0)
expect_equal(pd_cleaned, 
             pd_cleaned_frac[, 
                             colnames(pd_cleaned_frac) != "Fraction", 
                             with = FALSE])
# PD-TMT
pdtmt_input = data.table::fread("./raw_data/PDTMT/pdtmt_input.csv")
pdtmt_input2 = data.table::copy(pdtmt_input)
pdtmt_input2$ProteinMasterAccessions = pdtmt_input2$Protein.Accessions
pdtmt_input2$X..Protein.Groups = pdtmt_input2$X..Proteins
pdtmt_input3 = data.table::copy(pdtmt_input)
pdtmt_input3$Quan.Info = as.character(pdtmt_input3$Quan.Info)
pdtmt_input3$Quan.Info = "UNIQUE"
pdtmt_import = MSstatsConvert::MSstatsImport(list(input = pdtmt_input),
                                             "MSstatsTMT", "ProteomeDiscoverer")
pdtmt_import2 = MSstatsConvert::MSstatsImport(list(input = pdtmt_input2),
                                              "MSstatsTMT", "ProteomeDiscoverer")
pdtmt_import3 = MSstatsConvert::MSstatsImport(list(input = pdtmt_input3),
                                              "MSstatsTMT", "ProteomeDiscoverer")
pdtmt_cleaned = MSstatsConvert::MSstatsClean(pdtmt_import, 
                                             protein_id_column = "ProteinAccessions",
                                             remove_shared = TRUE)
pdtmt_cleaned2 = MSstatsConvert::MSstatsClean(pdtmt_import2, 
                                              protein_id_column = "ProteinMasterAccessions",
                                              remove_shared = TRUE)
pdtmt_cleaned3 = MSstatsConvert::MSstatsClean(pdtmt_import3, 
                                              protein_id_column = "ProteinMasterAccessions",
                                              remove_shared = TRUE)
expect_equal(pdtmt_cleaned, pdtmt_cleaned2)
expect_equal(pdtmt_cleaned[, colnames(pdtmt_cleaned) != "QuanInfo", 
                           with = FALSE], 
             pdtmt_cleaned3[, colnames(pdtmt_cleaned) != "QuanInfo", 
                            with = FALSE])
expect_error(MSstatsConvert::MSstatsClean(pdtmt_import, 
                                          protein_id_column = "ProteinAccessions",
                                          remove_shared = TRUE,
                                          intensity_columns_regexp = "Nothing"))
expect_equal(
    ncol(pdtmt_cleaned), 11
)
expect_true(nrow(pdtmt_cleaned) > 0)
# Progenesis
progenesis_input = data.table::fread("./raw_data/Progenesis/progenesis_input.csv")
progenesis_import = MSstatsConvert::MSstatsImport(list(input = progenesis_input),
                                                  "MSstats", "Progenesis")
progenesis_import2 = MSstatsConvert::MSstatsImport(
    list(input = rbind(progenesis_input[1, ], progenesis_input)),
    "MSstats", "Progenesis")
runs = unique(data.table::fread("./raw_data/Progenesis/progenesis_annot.csv")$Run)
pg_cleaned = MSstatsConvert::MSstatsClean(progenesis_import, runs)
pg_cleaned2 = MSstatsConvert::MSstatsClean(progenesis_import2, runs) 
expect_equal(ncol(pg_cleaned), 5)
expect_true(nrow(pg_cleaned) > 0)
expect_equal(pg_cleaned, pg_cleaned2)
# Skyline
skyline_input = data.table::fread("./raw_data/Skyline/skyline_input.csv")
skyline_input2 = data.table::copy(skyline_input)
skyline_input2$DetectionQValue = 0.01
skyline_input2$Truncated = ifelse(skyline_input2$Truncated, "True", "False")
skyline_input3 = data.table::copy(skyline_input)
skyline_input3$Fragment.Ion = rep(c("precursor", "a", "b"), times = nrow(skyline_input) / 3)
skyline_import = MSstatsConvert::MSstatsImport(list(input = skyline_input), 
                                               "MSstats", "Skyline")
skyline_import2 = MSstatsConvert::MSstatsImport(list(input = skyline_input2), 
                                                "MSstats", "Skyline")
skyline_import3 = MSstatsConvert::MSstatsImport(list(input = skyline_input3), 
                                                "MSstats", "Skyline")
sl_cleaned = MSstatsConvert::MSstatsClean(skyline_import)
sl_cleaned2 = MSstatsConvert::MSstatsClean(skyline_import2)
expect_equal(ncol(sl_cleaned), 13)
expect_true(nrow(sl_cleaned) > 0)
expect_equal(sl_cleaned, sl_cleaned2[, colnames(sl_cleaned), with = FALSE])
expect_error(MSstatsConvert:::.checkDDA(MSstatsConvert::MSstatsClean(skyline_import3)))
# SpectroMine
spectromine_input = data.table::fread("./raw_data/SpectroMine/spectromine_input.csv")
spectromine_import = MSstatsConvert::MSstatsImport(list(input = spectromine_input), 
                                                   "MSstatsTMT", "SpectroMine")
spectromine_import_error = MSstatsConvert::MSstatsImport(list(input = spectromine_input[, 1:19]),
                                                         "MSstatsTMT", "SpectroMine")
sm_cleaned = MSstatsConvert::MSstatsClean(spectromine_import)
expect_equal(ncol(sm_cleaned), 9)
expect_true(nrow(sm_cleaned) > 0)
expect_error(MSstatsConvert::MSstatsClean(spectromine_import_error))
# Spectronaut
spectronaut_input = data.table::fread("./raw_data/Spectronaut/spectronaut_input.csv")
spectronaut_input2 = data.table::copy(spectronaut_input)
spectronaut_input2$F.ExcludedFromQuantification = ifelse(
    spectronaut_input2$F.ExcludedFromQuantification,
    "True", "False"
)
spectronaut_import = MSstatsConvert::MSstatsImport(list(input = spectronaut_input), 
                                                   "MSstats", "Spectronaut")
spectronaut_import2 = MSstatsConvert::MSstatsImport(list(input = spectronaut_input2), 
                                                    "MSstats", "Spectronaut")
sn_cleaned = MSstatsConvert::MSstatsClean(spectronaut_import,
                                          intensity = "PeakArea")
sn_cleaned2 = MSstatsConvert::MSstatsClean(spectronaut_import2,
                                           intensity = "PeakArea")
expect_equal(ncol(sn_cleaned), 11)
expect_true(nrow(sn_cleaned) > 0)
expect_equal(sn_cleaned, sn_cleaned2)

# Metamorpheus
metamorpheus_table = data.table::fread("./raw_data/Metamorpheus/AllQuantifiedPeaks.tsv")
input = MSstatsConvert::MSstatsImport(list(input = metamorpheus_table), "MSstats", "Metamorpheus")
expect_identical(is(input), c("MSstatsMetamorpheusFiles", "MSstatsInputFiles"))