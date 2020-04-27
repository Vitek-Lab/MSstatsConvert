# TODO: ALL DATA CHECKS!!
.makeInclusionErrorMessage = function(parameter_name, legal_values, 
                                      information) {
    if(is.null(information)) {
        paste("Parameter", parameter_name, "must be one of", 
              paste(legal_values, sep = ", ", collapse = ", "))
    } else {
        paste(information, paste(legal_values, sep = ", ", collapse = ", "))
    }
}

.isLegalValue = function(parameter, legal_values = NULL, 
                         can_be_null = FALSE, information = NULL) {
    parameter_name = deparse(substitute(parameter))
    if(is.null(parameter)) {
        if(!can_be_null) {
            stop(paste("Parameter", parameter_name, "cannot be NULL"))    
        }
    } else {
        if(!is.null(legal_values)) {
            if(!is.element(parameter, legal_values)) {
                stop(.makeInclusionErrorMessage(parameter_name, legal_values,
                                                information))
            }    
        }
    }
    parameter
}


.checkConverterParams = function(remove_shared,
                                 remove_few_measurements,
                                 summary_for_multiple_rows, 
                                 handle_single_feature_per_protein, ...) {
    inputs = list(...)
    if (!all(sapply(inputs, function(x) inherits(x, "data.frame") | is.character(x)))) {
        msg = "Inputs and annotation must be either a data.frame / data.table or a path to a file"
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    if (!(identical(summary_for_multiple_rows, max) | 
        identical(summary_for_multiple_rows, sum))) {
        msg = "Summary for multiple rows must be one of sum, max"
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    if (!is.logical(c(remove_shared, remove_few_measurements))) {
        msg = "Remove shared, ... , has to be logical (TRUE or FALSE)"
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    }
    
}


# checks for data? 
# could have three functions: .checkInput, .checkAnnotation, .checkParams


#' Select an available options from a set of possibilities
#' @param possibilities possible legal values of a variable
#' @param option_set set of values that includes one of the `possibilities`
#' @param fall_back if there is none of the `possibilities` in `option_set`,
#' default to `fall_back`
#' @return same as option_set, usually character
.findAvailable = function(possibilities, option_set, fall_back = NULL) {
    chosen = option_set[option_set %in% possibilities]
    if (length(chosen) != 1L) {
        if (is.null(fall_back)) {
            NULL
        } else {
            if (fall_back %in% possibilities) {
                fall_back
            } else {
                NULL
            }
        }
    } else {
        chosen
    }
    # TODO: throw + log error here if couldn't find
    # TODO: there should be a specific value that will be chosen if it's 
    # among chosen
}
