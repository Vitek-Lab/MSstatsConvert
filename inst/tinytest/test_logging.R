# Null appender
expect_null(
    MSstatsConvert:::.nullAppender("INFO", "MSG")
)
# For .onLoad
expect_true(!is.null(getOption("MSstatsLog")))
expect_true(!is.null(getOption("MSstatsMsg")))
