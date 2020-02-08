old_1 = c("Protein", "Peptide")
old_2 = c("Protein", "Peptide", "Other")
update = c("Protein" = "ProteinName", 
           "Peptide" = "PeptideSequence")
x <- data.frame('Protein' = 'x', 'Peptide' = 'y', 'Other' = 'Other')

test_that("Column update works", {
    expect_equal(old_1, 
                 .updateColnames(x[, 1:2], 
                                 c("Protein" = "Protein", "Peptide" = "Peptide")))
    expect_equal(c("ProteinName", "PeptideSequence"),
                 .updateColnames(x[, 1:2], update))
    expect_equal(c("ProteinName", "PeptideSequence", "Other"),
                 .updateColnames(x, update))
})

