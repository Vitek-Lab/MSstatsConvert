# Null appender
expect_null(
    MSstatsConvert:::.nullAppender("INFO", "MSG")
)
# For .onLoad
expect_true(!is.null(getOption("MSstatsLog")))
expect_true(!is.null(getOption("MSstatsMsg")))

old_opts = options("MSstatsLog", "MSstatsMsg")

f = tempfile(fileext = ".log")
writeLines("old line", f)
MSstatsConvert::MSstatsLogsSettings(TRUE, append = FALSE, verbose = FALSE,
                                    log_file_path = f)
getOption("MSstatsLog")("INFO", "first")
getOption("MSstatsLog")("INFO", "second")
lines = readLines(f)
expect_equal(length(lines), 2L)
expect_false(any(grepl("old line", lines)))

MSstatsConvert::MSstatsLogsSettings(TRUE, append = TRUE, verbose = FALSE,
                                    log_file_path = f)
getOption("MSstatsLog")("INFO", "third")
expect_equal(length(readLines(f)), 3L)

options(old_opts)
