#' log4r appender used not to write messages 
#' 
#' A convenience function written to save time on checking if messages should
#' be printed or logs should be written to a file.
#' 
#' @param level log level 
#' @param ... messages - ignored
#' @keywords internal
.nullAppender = function(level, ...) {
    invisible(NULL)
}


#' Set default logging object when package is loaded
#' @param ... ignored
#' @importFrom log4r file_appender console_appender
#' @keywords internal
.onLoad = function(...) {
    logs = getOption("MSstatsLog")
    msgs = getOption("MSstatsMsg")
    time_now = Sys.time()
    path = paste0("./MSstats_log_", gsub("[ :\\-]", "_", time_now), ".log")
    
    if (is.null(logs)) {
        ms_logs = file_appender(path, append = TRUE)
        options(MSstatsLog = ms_logs)
    }
    if (is.null(msgs)) {
        ms_messages = console_appender()
        options(MSstatsMsg = ms_messages)
    }
    
    if (is.null(logs)) {
        getOption("MSstatsLog")("INFO", paste("Initialized MSstats session:",
                                              time_now))
    }
}


#' Set log4 appenders based on user preferences
#' @param log logical, if TRUE, messages will be written to a file.
#' @param append logical, if TRUE, messages will be written to an existing 
#' MSstats log file.
#' @param verbose logical, if TRUE, messages will be printed to the console.
#' @importFrom log4r file_appender
#' @keywords internal
.setMSstatsLogger = function(log, append, verbose) {
    if (log) {
        if (append) {
            msstats_logs = sort(list.files(".", "MSstats_log"),
                                decreasing = TRUE)
            if(length(msstats_logs) > 0) {
                options(MSstatsLog = file_appender(msstats_logs[1], 
                                                   append = TRUE))
            }
        }
    } else {
        options(MSstatsLog = .nullAppender)
    }
    
    if (!verbose) {
        options(MSstatsMsg = .nullAppender)
    }
}
