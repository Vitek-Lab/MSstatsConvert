# Test intensity parameter
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$FG.MS1Quantity = 100000
msstats_input = MSstatsConvert::MSstatsImport(
    list(input = spectronaut_raw), "MSstats", "Spectronaut")
output = MSstatsConvert:::.cleanRawSpectronaut(msstats_input, intensity = 'FG.MS1Quantity', 
                                    calculateAnomalyScores = FALSE, 
                                    anomalyModelFeatures = c())
expect_true(all(output$Intensity == 100000))
expect_false("Fraction" %in% colnames(output))

spectronaut_frac = data.table::copy(spectronaut_raw)
spectronaut_frac$R.Fraction = 2L
msstats_frac = MSstatsConvert::MSstatsImport(
    list(input = spectronaut_frac), "MSstats", "Spectronaut")
output_frac = MSstatsConvert:::.cleanRawSpectronaut(msstats_frac, intensity = 'FG.MS1Quantity',
                                    calculateAnomalyScores = FALSE,
                                    anomalyModelFeatures = c())
expect_true("Fraction" %in% colnames(output_frac))
expect_false("RFraction" %in% colnames(output_frac))
expect_true(all(output_frac$Fraction == 2L))

spectronaut_frac_na = data.table::copy(spectronaut_raw)
spectronaut_frac_na$R.Fraction = NA_integer_
msstats_frac_na = MSstatsConvert::MSstatsImport(
    list(input = spectronaut_frac_na), "MSstats", "Spectronaut")
output_frac_na = MSstatsConvert:::.cleanRawSpectronaut(msstats_frac_na, intensity = 'FG.MS1Quantity',
                                    calculateAnomalyScores = FALSE,
                                    anomalyModelFeatures = c())
expect_true("Fraction" %in% colnames(output_frac_na))
expect_false(any(is.na(output_frac_na$Fraction)))
expect_true(all(output_frac_na$Fraction == 1))

expect_error(MSstatsConvert:::.cleanRawSpectronaut(msstats_input, intensity = 'invalid',
                                                   calculateAnomalyScores = FALSE,
                                                   anomalyModelFeatures = c()),
             pattern = "not found in input data")

# Tests for .assignSpectronautIsotopeLabelType
make_spec_input = function(peptides) {
    data.table::data.table(PeptideSequence = peptides)
}

dt = make_spec_input(c(
    "_PEPTIDEK[Lys6]_",   # heavy via Lys6
    "_PEPTIDER[Arg10]_",  # heavy via Arg10
    "_PEPTIDEK_",         # light (has K, but no Lys6 label)
    "_PEPTIDER_",         # light (has R, but no Arg10 label)
    "_PEPTIDEAC_"         # NA (no K or R)
))
result = MSstatsConvert:::.assignSpectronautIsotopeLabelType(
    dt, heavyLabels = c("K[Lys6]", "R[Arg10]"))
expect_equal(result$IsotopeLabelType, c("H", "H", "L", "L", NA_character_))
expect_equal(result$PeptideSequence,
             c("_PEPTIDEK_", "_PEPTIDER_", "_PEPTIDEK_", "_PEPTIDER_", "_PEPTIDEAC_"))

dt = make_spec_input(c(
    "_PEPTM[Oxidation]IDEK_",   
    "_S[Kmodification]PEPTIDEK_", 
    "_PEPTM[Oxidation]IDEK[Lys6]_",
    "_ACDEGFHI_"
))
result = MSstatsConvert:::.assignSpectronautIsotopeLabelType(
    dt, heavyLabels = "K[Lys6]")
expect_equal(result$IsotopeLabelType, c("L", "L", "H", NA_character_))
# Non-heavy brackets and heavy label brackets are stripped from PeptideSequence
expect_equal(result$PeptideSequence,
             c("_PEPTM[Oxidation]IDEK_", "_S[Kmodification]PEPTIDEK_",
               "_PEPTM[Oxidation]IDEK_", "_ACDEGFHI_"))

dt = make_spec_input(c("_PEPTIDEK_", "_PEPTIDER_"))
result = MSstatsConvert:::.assignSpectronautIsotopeLabelType(dt, heavyLabels = NULL)
expect_equal(result, dt)
