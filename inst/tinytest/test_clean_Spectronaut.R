# Test intensity parameter
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$FG.MS1Quantity = 100000
msstats_input = MSstatsConvert::MSstatsImport(
    list(input = spectronaut_raw), "MSstats", "Spectronaut")
output = MSstatsConvert:::.cleanRawSpectronaut(msstats_input, intensity = 'MS1Quantity', 
                                    calculateAnomalyScores = FALSE, 
                                    anomalyModelFeatures = c())
expect_true(all(output$Intensity == 100000))

expect_error(MSstatsConvert:::.cleanRawSpectronaut(msstats_input, intensity = 'invalid', 
                                                   calculateAnomalyScores = FALSE, 
                                                   anomalyModelFeatures = c()), 
             pattern = "'arg' should be one of .*PeakArea")
