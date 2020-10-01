openms = system.file("tinytest/raw_data/OpenMS/openms_input.csv",
                     package = "MSstatsConvert")
openms = read.csv(openms)

imported = MSstatsConvert::MSstatsImport(list(input = openms),
                                         "MSstats", "OpenMS")
tinytest::expect_identical(is(imported), 
                           c("MSstatsOpenMSFiles", "MSstatsInputFiles"))
cleaned = MSstatsConvert::MSstatsClean(imported)
tinytest::expect_identical(class(cleaned), 
                           c("data.table", "data.frame"))
annotation = MSstatsConvert::MSstatsMakeAnnotation(cleaned, NULL)
tinytest::expect_identical(class(annotation), 
                           c("data.table", "data.frame"))
processed = MSstatsConvert::MSstatsPreprocess(
    cleaned, annotation, 
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
)
tinytest::expect_identical(class(processed), 
                           c("data.table", "data.frame"))
tinytest::expect_true(nrow(processed) > 0)
tinytest::expect_true(ncol(processed) == 10)
balanced = MSstatsConvert::MSstatsBalancedDesign(
    processed, 
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
)
tinytest::expect_true(nrow(balanced) > 0)
tinytest::expect_true(ncol(balanced) == 11)
