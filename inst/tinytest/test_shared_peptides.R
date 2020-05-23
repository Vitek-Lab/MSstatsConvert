data_with_shared = data.table::data.table(
    ProteinName = c("A", "A", "B", "B", "B", 
                    "C", "D", "D", "D", "D"),
    PeptideSequence = c("P1", "P2", "P1", "P3", "P4",
                        "P5", "P6", "P3", "P7", "P8"))
data_no_shared = data.table::data.table(ProteinName = letters[1:10],
                                        PeptideSequence = letters[1:10])

remove_shared = MSstatsConvert:::.handleSharedPeptides(data_with_shared, TRUE)
dont_remove_shared = MSstatsConvert:::.handleSharedPeptides(data_with_shared, FALSE)

no_shared = MSstatsConvert:::.handleSharedPeptides(data_no_shared)

# Correct output class
expect_equal(class(remove_shared)[1], "data.table")
# Nothing changes if remove = FALSE
expect_identical(data_with_shared, dont_remove_shared)
# Nothing changes if there are no shared peptides
expect_identical(no_shared, data_no_shared)
# Shared peptides are removed
expect_equal(nrow(remove_shared), 6L)
expect_true(all(!(c("P1", "P3") %in% remove_shared$PeptideSequence)))
