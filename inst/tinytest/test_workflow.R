openms = system.file("tinytest/raw_data/OpenMS/openms_input.csv",
                     package = "MSstatsConvert")
openms = read.csv(openms)

imported = MSstatsConvert::MSstatsImport(list(input = openms),
                                         "MSstats", "OpenMS")
expect_identical(is(imported), 
                           c("MSstatsOpenMSFiles", "MSstatsInputFiles"))
cleaned = MSstatsConvert::MSstatsClean(imported)
expect_identical(class(cleaned), 
                           c("data.table", "data.frame"))
annotation = MSstatsConvert::MSstatsMakeAnnotation(cleaned, NULL)
expect_identical(class(annotation), 
                           c("data.table", "data.frame"))
processed = MSstatsConvert::MSstatsPreprocess(
    cleaned, annotation, 
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
)
expect_identical(class(processed), 
                           c("data.table", "data.frame"))
expect_true(nrow(processed) > 0)
expect_true(ncol(processed) == 10)
balanced = MSstatsConvert::MSstatsBalancedDesign(
    processed, 
    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
)
expect_true(nrow(balanced) > 0)
expect_true(ncol(balanced) == 11)
# Placeholder with simple test
expect_equal(1 + 1, 3)
