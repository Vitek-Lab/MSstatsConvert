#' log4r appender used not to write messages 
#' 
#' A convenience function written to save time on checking if messages should
#' be printed or logs should be written to a file.
#' 
#' @param level log level 
#' @param ... messages - ignored
#' @return NULL invisibly
#' @keywords internal
.nullAppender = function(level, ...) {
    invisible(NULL)
}


#' Set default logging object when package is loaded
#' @param ... ignored
#' @importFrom log4r file_appender console_appender
#' @return none, sets options called MSstatsLog and MSstatsMsg
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
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' @param base start of the file name.
#' 
#' @return TRUE invisibly in case of successful logging setup.
#' @export
#' 
#' @examples 
#' # No logging and no messages
#' MSstatsLogsSettings(FALSE, FALSE, FALSE)
#' # Log, but do not display messages
#' MSstatsLogsSettings(TRUE, FALSE, FALSE)
#' # Log to an existing file
#' file.create("new_log.log")
#' MSstatsLogsSettings(TRUE, TRUE, log_file_path = "new_log.log")
#' # Do not log, but display messages
#' MSstatsLogsSettings(FALSE)
#' 
MSstatsLogsSettings = function(use_log_file = TRUE, append = FALSE,
                               verbose = TRUE, log_file_path = NULL,
                               base = "MSstats_log_") {
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
                path = paste0(base, gsub("[ :\\-]", "_", time_now), 
                              ".log")
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


#' Make a message about filtering based on fixed values
#' @inheritParams .filterExact
#' @return character - message
#' @keywords internal
.makeExactFilterMessage = function(col_name, filter_symbols, 
                                   behavior, fill_value
) {
   if (behavior == "remove") {
       subject = "Rows"
       event = "are removed"
       what = ""
   } else {
       subject = "Intensities"
       event = "are replaced with"
       what = fill_value
   }
    paste("**", subject, "with values of", col_name, "equal to", 
          paste(filter_symbols, sep = ", ", collapse = ", "),
          event, what)
}


#' Make a message about filtering based on a score
#' @inheritParams .filterByScore
#' @return character - message
#' @keywords internal
.makeScoreFilterMessage = function(score_column, score_threshold, direction,
                                   behavior, fill_value) {
    if (behavior == "remove") {
        subject = "Rows"
        event = "are removed"
        what = ""
    } else {
        subject = "Intensities"
        event = "are replaced with"
        what = fill_value
    } 
    msg = paste("**", subject, "with values", direction, "than", score_threshold,
                "in", score_column, event, what)
    msg
}


#' Log information about converter options
#' @inheritParams MSstatsPreprocess
#' @return TRUE invisibly if message was logged 
#' @keywords internal
.logConverterOptions = function(feature_columns, remove_shared_peptides,
                                        remove_single_feature_proteins,
                                        feature_cleaning) {
    init = "** The following options are used:"
    features = paste("  - Features will be defined by the columns:",
                     paste(feature_columns, sep = ", ", collapse = ", "))
    if (remove_shared_peptides) {
        shared = "  - Shared peptides will be removed"
    } else {
        shared = "  - Shared peptides will not be removed"
    }
    if (remove_single_feature_proteins) {
        single = "  - Proteins with a single feature will be removed"
    } else {
        single = "  - Proteins with single feature will not be removed"
    }
    if (feature_cleaning[["remove_features_with_few_measurements"]]) {
        few = "  - Features with less than 3 measurements will be removed"
    } else {
        few = "  - Features with less than 3 measurements will be kept"
    }
    msg = paste(init, features, shared, single, few, sep = "\n")
    getOption("MSstatsLog")("INFO", msg)    
    getOption("MSstatsMsg")("INFO", msg)
    invisible(TRUE)
}


#' Make a message about successful data cleaning/importing
#' @param tool name of a signal processing tool
#' @return TRUE invisibly if logging was sucessful
#' @keywords internal
.logSuccess = function(tool, event) {
    if (event == "clean") {
        what = "cleaned succesfully."
    } else {
        what = "imported succesfully."
    }
    msg = paste("** Raw data from", tool, what)
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    invisible(TRUE)
}

#' Save session information
#' 
#' @param path optional path to output file. If not provided, "MSstats_session_info" 
#' and current timestamp will be used as a file name
#' @param append if TRUE and file given by the `path` parameter already exists,
#' session info will be appended to the file
#' @return TRUE invisibly after session info was saved
#' @export
#' @importFrom utils sessionInfo
#' @examples 
#' MSstatsSaveSessionInfo("session_info.txt")
#' 
MSstatsSaveSessionInfo = function(path = NULL, append = TRUE) {
    if (is.null(path)) {
        time_now = Sys.time()
        path = paste0("./MSstats_session_info_", 
                      gsub("[ :\\-]", "_", time_now), ".log")
    }
    session_info = utils::sessionInfo()
    sink(path, append = append)
    print(session_info)
    sink()
    invisible(TRUE)
}
