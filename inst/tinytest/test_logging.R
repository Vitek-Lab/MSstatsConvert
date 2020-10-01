# Null appender
tinytest::expect_null(
    MSstatsConvert:::.nullAppender("INFO", "MSG")
)
# For .onLoad
tinytest::expect_true(!is.null(getOption("MSstatsLog")))
tinytest::expect_true(!is.null(getOption("MSstatsMsg")))
