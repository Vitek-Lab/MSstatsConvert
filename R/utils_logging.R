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
        ms_logs = file_appender(path)
        options(MSstatsLog = ms_logs)
    }
    if (is.null(msgs)) {
        ms_messages = console_appender()
        options(MSstatsMsg = ms_messages)
    }
}


#' Set how MSstats will log information from data processing
#' 
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be added
#' to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be printed
#' to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' 
#' @return TRUE invisibly in case of successful logging setup.
#' @export
#' 
MSstatsLogsSettings = function(use_log_file = TRUE, append = FALSE,
                               verbose = TRUE, log_file_path = NULL) {
    checkmate::assertLogical(use_log_file)
    checkmate::assertLogical(append)
    checkmate::assertLogical(verbose)
    checkmate::assertCharacter(log_file_path, null.ok = TRUE)
    checkmate::assertTRUE(!(append & !use_log_file))
    checkmate::assertTRUE(!(append & is.null(log_file_path)))
    if (!is.null(log_file_path)) {
        checkmate::assertPathForOutput(log_file_path, overwrite = TRUE)
    }
    
    if (use_log_file) {
        if (append) {
            options("MSstatsLog" = log4r::file_appender(log_file_path))
        } else {
            if (!is.null(log_file_path)) {
                options("MSstatsLog" = log4r::file_appender(log_file_path),
                        append = append)
            } else {
                time_now = Sys.time()
                path = paste0("./MSstats_log_", gsub("[ :\\-]", "_", time_now), ".log")
                options("MSstatsLog" = log4r::file_appender(path))
            }
        }
    } else {
        options("MSstatsLog" = .nullAppender)
    }
    
    if (verbose) {
        options("MSstatsMsg" = console_appender())
    } else {
        options("MSstatsMsg" = .nullAppender)
    }
    invisible(TRUE)
}
