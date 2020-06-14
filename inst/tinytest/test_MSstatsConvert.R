data_invalid = data.table::data.table(ProteinName = letters[1:10],
                                         PeptideSequence = LETTERS[1:10])
data_no_channel = data.table::data.table(ProteinName = letters[1:10],
                                         PeptideSequence = LETTERS[1:10],
                                         IsotopeLabelType = "L")
data_no_channel_LH = data.table::data.table(ProteinName = letters[1:10],
                                            PeptideSequence = LETTERS[1:10],
                                            IsotopeLabelType = rep(c("L", "H"),
                                                                   each = 5))
data_channel = data.table::data.table(ProteinName = letters[1:10],
                                      PeptideSequence = LETTERS[1:10],
                                      Channel = paste0("A", 1:10))
tmt_object = MSstatsConvert:::.MSstatsFormat(data_channel)
labelfree_object = MSstatsConvert:::.MSstatsFormat(data_no_channel)
labelled_object = MSstatsConvert:::.MSstatsFormat(data_no_channel_LH)

# Has the right format for MSstats data
expect_true(is(labelled_object, "MSstatsValidated"))
expect_true(is(labelfree_object, "MSstatsValidated"))
expect_true(is(tmt_object, "MSstatsValidated"))
# Inherits after data.frame
expect_true(is(labelled_object, "data.frame"))
expect_true(is(labelfree_object, "data.frame"))
expect_true(is(tmt_object, "data.frame"))