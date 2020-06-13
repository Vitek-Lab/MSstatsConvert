test_path = system.file("tinytest", "get_dt_test.csv", package = "MSstatsConvert")
test_df = data.table::fread(test_path)

# Import and conversion ----
## data.table from path
tinytest::expect_equal(class(MSstatsConvert:::.getDataTable(test_path))[1], "data.table")
## data.table from data.frame
tinytest::expect_equal(class(MSstatsConvert:::.getDataTable(test_df))[1], "data.table")
## When it's not a data.frame, a test 
tinytest::expect_error(MSstatsConvert:::.getDataTable(numeric(10)))
# Column names ----
## Update colnames when needed
tinytest::expect_equal(MSstatsConvert:::.updateColnames(test_df, "x", "z"), c("z", "y"))
## Do not update, when the column doesn't exist
tinytest::expect_equal(MSstatsConvert:::.updateColnames(test_df, "z", "z"), c("x", "y"))
## Remove unnecessary symbols from column names
tinytest::expect_equal(MSstatsConvert:::.standardizeColnames("Col.name"), "Colname")
tinytest::expect_equal(MSstatsConvert:::.standardizeColnames("[Col%name]"), "Colname")
### Does not cover ":" - is it bad?
tinytest::expect_equal(MSstatsConvert:::.getChannelColumns(c("Ab:1", "Ab:2"), "Ab"), c("Ab:1", "Ab:2"))
# Additional columns
updated_df = MSstatsConvert:::.fillValues(test_df, list(z = NA))
tinytest::expect_true(all(is.na(updated_df$z)) & is.element("z", colnames(updated_df)))
# Column types ---
updated_col_types = updated_df[, list(x = as.character(x), y = as.factor(y), z = as.numeric(z))]
fixed_column_types = MSstatsConvert:::.fixColumnTypes(updated_col_types, numeric_columns = "x", 
                                                      character_columns = "y",
                                                      factor_columns = "z")
tinytest::expect_equal(class(fixed_column_types$x), "numeric")
tinytest::expect_equal(class(fixed_column_types$y), "character")
tinytest::expect_equal(class(fixed_column_types$z), "factor")
