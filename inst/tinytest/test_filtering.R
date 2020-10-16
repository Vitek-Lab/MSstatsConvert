# Ceate data ----
data_to_filter = data.table::data.table(
    ProteinName = rep(LETTERS[1:5], 10),
    PeptideSequence = rep(letters[1:10], 5),
    PeptideSequenceModified = rep(letters[1:10], 5),
    Intensity = seq(0, 1, length.out = 50),
    Score = seq(0, 1, length.out = 50),
    Symbol_1 = c(rep("+", 25), rep("-", 25)),
    Symbol_2 = c(rep("+", 25), rep("-", 25))
)
data_to_filter$Score[c(1, 13, 15, 27, 38)] = NA
# Filter by pattern ----
## Warn when the column isn't there
expect_stdout(
    MSstatsConvert:::.filterByPattern(data_to_filter, "Not_there", "\\+", TRUE, FALSE)
)
## Filter and don't drop the column
expect_identical(
    MSstatsConvert:::.filterByPattern(data_to_filter, "Symbol_1", "\\+", TRUE, FALSE),
    data_to_filter[26:50, ]
)
## Filter and drop the column
expect_identical(
    MSstatsConvert:::.filterByPattern(data_to_filter, "Symbol_1", "\\+", TRUE, TRUE),
    data_to_filter[26:50, -6]
)
## Don't filter and don't drop the column
expect_identical(
    MSstatsConvert:::.filterByPattern(data_to_filter, "Symbol_1", "\\+", FALSE, FALSE),
    data_to_filter
)
# Don't filter but drop the column 
expect_identical(
    MSstatsConvert:::.filterByPattern(data_to_filter, "Symbol_1", "\\+", FALSE, TRUE),
    data_to_filter[, -6]
)
# Filter by exact values ----
## Warn when the column isn't there
expect_stdout(MSstatsConvert:::.filterExact(data_to_filter, "Not_there", 
                                             "-", "remove", NULL, TRUE, TRUE))
## Filter and drop column
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               "-", "remove", NULL, TRUE, TRUE),
                 data_to_filter[1:25, -6])
## Filter so that nothing's left
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               c("+", "-"), "remove", NULL, TRUE, TRUE),
                 data_to_filter[FALSE, -6])
## Filter, but don't drop the column
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               "+", "remove", NULL, TRUE, FALSE),
                 data_to_filter[26:50, ])
## Don't filter, but drop column
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               c("+", "-"), "remove", NULL, FALSE, TRUE),
                 data_to_filter[, -6])
## Don't filter and don't drop column
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               c("+", "-"), "remove", NULL, FALSE, FALSE),
                 data_to_filter)
## Filter and drop, but the symbol is not there
expect_identical(MSstatsConvert:::.filterExact(data_to_filter, "Symbol_1", 
                                               "R", "remove", 
                                               NULL, TRUE, TRUE),
                 data_to_filter[, -6])
# Filter multiple columns by exact values ----
## Warn when the column isn't there
expect_stdout(
    MSstatsConvert:::.filterManyColumns(data_to_filter, c("Not_there", "Not_there"), "+")
)
## Filter
expect_identical(
    MSstatsConvert:::.filterManyColumns(data_to_filter, c("Symbol_1", "Symbol_2"), "+"),
    data_to_filter[26:50, 1:5]
)
## When the values aren't there just drop columns
expect_identical(
    MSstatsConvert:::.filterManyColumns(data_to_filter, c("Symbol_1", "Symbol_2"), "X"),
    data_to_filter[, 1:5]
)
# Filter by a numerical score ----
## Warn when the column is not there
expect_stdout(MSstatsConvert:::.filterByScore(data_to_filter, "Not_there", 0.3, "smaller",
                                               "fill", handle_na = "keep", fill_value = NA,
                                               filter = TRUE, drop = TRUE)
)
## Don't filter and don't drop
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = FALSE, drop = FALSE),
                 data_to_filter)
## Don't filter but drop
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = FALSE, drop = TRUE),
                 data_to_filter[, -5])
## Keep only scores bigger than 0.3, don't drop the column, keep NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = TRUE, drop = FALSE),
                 data_to_filter[Score > 0.3 | is.na(Score), ])
## Keep only scores bigger than 0.3, don't drop the column, remove NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "remove", fill_value = NA,
                                                 filter = TRUE, drop = FALSE),
                 data_to_filter[Score > 0.3 & !is.na(Score), ])
## Keep only scores bigger than 0.3, drop the column, keep NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = TRUE, drop = TRUE),
                 data_to_filter[Score > 0.3 | is.na(Score), -5])
## Keep only scores bigger than 0.3, drop the column, remove NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                 "remove", handle_na = "remove", fill_value = NA,
                                                 filter = TRUE, drop = TRUE),
                 data_to_filter[Score > 0.3 & !is.na(Score), -5])
## Keep only scores smaller than 0.3, don't drop the column, keep NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = TRUE, drop = FALSE),
                 data_to_filter[Score < 0.3 | is.na(Score), ])
## Keep only scores smaller than 0.3, don't drop the column, remove NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                 "remove", handle_na = "remove", fill_value = NA,
                                                 filter = TRUE, drop = FALSE),
                 data_to_filter[Score < 0.3 & !is.na(Score), ])
## Keep only scores smaller than 0.3, drop the column, keep NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                 "remove", handle_na = "keep", fill_value = NA,
                                                 filter = TRUE, drop = TRUE),
                 data_to_filter[Score < 0.3 | is.na(Score), -5])
## Keep only scores smaller than 0.3, drop the column, remove NA
expect_identical(MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                 "remove", handle_na = "remove", fill_value = NA,
                                                 filter = TRUE, drop = TRUE),
                 data_to_filter[Score < 0.3 & !is.na(Score), -5])
## Keep only scores greater than 0.3, replace other Intensities with NA,
## drop the column, keep NA
fill_smaller_than_0_3 = MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                        "fill", handle_na = "keep", fill_value = NA,
                                                        filter = TRUE, drop = TRUE)
expect_identical(ncol(fill_smaller_than_0_3), 6L)
expect_identical(fill_smaller_than_0_3$Intensity, 
                 ifelse(data_to_filter$Score > 0.3 | is.na(data_to_filter$Score),
                        data_to_filter$Intensity, NA))
## Keep only scores bigger than 0.3, drop the column, remove NA, 
## replace other Intensities with NA
fill_smaller_than_0_3_2 = MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "greater",
                                                          "fill", handle_na = "remove", fill_value = NA,
                                                          filter = TRUE, drop = TRUE)
expect_identical(ncol(fill_smaller_than_0_3_2), 6L)
expect_identical(fill_smaller_than_0_3_2$Intensity, 
                 ifelse(data_to_filter$Score > 0.3 & !is.na(data_to_filter$Score),
                        data_to_filter$Intensity, NA))

## Keep only scores smaller than 0.3, replace other Intensities with NA,
## drop the column, keep NA
fill_smaller_than_0_3_s = MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                          "fill", handle_na = "keep", fill_value = NA,
                                                          filter = TRUE, drop = TRUE)
expect_identical(ncol(fill_smaller_than_0_3_s), 6L)
expect_identical(fill_smaller_than_0_3_s$Intensity, 
                 ifelse(data_to_filter$Score < 0.3 | is.na(data_to_filter$Score),
                        data_to_filter$Intensity, NA))
## Keep only scores smaller than 0.3, drop the column, remove NA, 
## replace other Intensities with NA
fill_smaller_than_0_3_2_s = MSstatsConvert:::.filterByScore(data_to_filter, "Score", 0.3, "smaller",
                                                            "fill", handle_na = "remove", fill_value = NA,
                                                            filter = TRUE, drop = TRUE)
expect_identical(ncol(fill_smaller_than_0_3_2_s), 6L)
expect_identical(fill_smaller_than_0_3_2_s$Intensity, 
                 ifelse(data_to_filter$Score < 0.3 & !is.na(data_to_filter$Score),
                        data_to_filter$Intensity, NA))
# Full filtering function ----
## Nothing to filter
expect_equal(
    MSstatsConvert:::.handleFiltering(data_to_filter, 
                                      list(), 
                                      list(), 
                                      list()),
    data_to_filter
)
## Filter and drop
expect_identical(
    MSstatsConvert:::.handleFiltering(data_to_filter, 
                                      list(list(score_column = "Score", score_threshold = 0.3, 
                                                direction = "greater", behavior = "remove", 
                                                handle_na = "keep", fill_value = NA,
                                                filter = FALSE, drop_column = TRUE)), 
                                      list(list(col_name = "Symbol_1", 
                                                filter_symbols = "+", filter = TRUE,
                                                behavior = "remove", fill_value = NULL,
                                                drop_column = TRUE)),
                                      list(list(col_name = "Symbol_2", pattern = "\\+", 
                                                filter = TRUE, drop_column = TRUE))),
    data_to_filter[(Score > 0.3 | is.na(Score)) & Symbol_1 == "-", 1:4]
)
## Filter, but don't drop
expect_equal(
    MSstatsConvert:::.handleFiltering(data_to_filter, 
                                      list(list(score_column = "Score", score_threshold = 0.3, 
                                                direction = "greater", behavior = "remove", 
                                                handle_na = "keep", fill_value = NA,
                                                filter = FALSE, drop_column = TRUE)), 
                                      list(list(col_name = "Symbol_1", 
                                                filter_symbols = "+", filter = TRUE, 
                                                behavior = "remove", fill_value = NULL,
                                                drop_column = FALSE)),
                                      list(list(col_name = "Symbol_1", pattern = "\\+", 
                                                filter = TRUE, drop_column = TRUE))),
    data_to_filter[(Score > 0.3 | is.na(Score)) & Symbol_1 == "-", c(1:4, 7)]
)
