#' Check if there are overlapping features and remove if needed
#' @param input data.table preprocessed by one of the .cleanRaw* functions and
#' merged with annotation.
#' @return `data.table`
#' @keywords internal
.handleFractions = function(input) {
    if (is.element("Channel", colnames(input))) {
        fractions = unique(input$Fraction)
        if (length(fractions) > 1) {
            input = .handleFractionsTMT(input)
            msg = "Fractions belonging to same mixture have been combined."
            getOption("MSstatsLog")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg)
        } else {
            input$Fraction = 1
        }
    } else {
        input = .handleFractionsLF(input)
    }
    input
}

#' Remove peptide ions overlapped among multiple fractions of the same biological mixture
#' @inheritParams .handleFractions
#' @return `data.table`
#' @keywords internal
.handleFractionsTMT = function(input) {
    techrun = feature = id = Mixture = TechRepMixture = Run = NULL
    
    input[, techrun := paste(Mixture, TechRepMixture, sep = "_")]
    input[, id := paste(feature, Run, sep = "_")]
    
    unoverlapped_list = vector("list", length(unique(input$techrun)))
    names(unoverlapped_list) = unique(input$techrun)
    for (technical_run in unique(input$techrun)) {
        single_run = input[techrun == technical_run & !is.na(Intensity), ]
        features_to_remove = .getOverlappingFeatures(single_run)
        if (length(features_to_remove) > 0) {
            single_run = .filterOverlapped(single_run, mean, features_to_remove)
            features_to_remove = .getOverlappingFeatures(single_run)
            msg = paste("For peptides overlapped between fractions of",
                        technical_run, "use the fraction with maximal average abundance.")
            getOption("MSstatsLog")("INFO", msg)
            getOption("MSstatsMsg")("INFO", msg)
            if (length(features_to_remove) > 0) {
                single_run = .filterOverlapped(single_run, sum, features_to_remove)
                features_to_remove = .getOverlappingFeatures(single_run)
                msg = paste("For peptides overlapped between fractions of",
                            technical_run, "use the fraction with maximal summation abundance.")
                getOption("MSstatsLog")("INFO", msg)
                getOption("MSstatsMsg")("INFO", msg)
                if (length(features_to_remove) > 0) {
                    single_run = .filterOverlapped(single_run, max, features_to_remove)
                    msg = paste("For peptides overlapped between fractions of",
                                technical_run, "use the fraction with maximal abundance.")
                    getOption("MSstatsLog")("INFO", msg)
                    getOption("MSstatsMsg")("INFO", msg)
                    if (data.table::uniqueN(input$Run) > 1) {
                        single_run = single_run[, list(Intensity = mean(Intensity, na.rm = TRUE)),
                                                by = setdiff(colnames(single_run), 
                                                             c("Run", "Intensity",
                                                               "Fraction", "id", "n_psms",
                                                               "QuanInfo", "IonsScore",
                                                               "IsolationInterference"))]
                    }
                }
            }
        }
        unoverlapped_list[[technical_run]] = single_run[, !(colnames(single_run) %in% c("run", "techrun", "id")),
                                                        with = FALSE]
    }
    input = rbindlist(unoverlapped_list)
    input[, Run := paste(Mixture, TechRepMixture, sep = "_")]
    input
}


#' Get features that are overlapped among multiple runs
#' @inheritParams .handleFractions
#' @return `data.table`
#' @keywords internal
.getOverlappingFeatures = function(input) {
    Run = feature = n_runs = NULL
    
    count_fractions = input[, list(n_runs = uniqueN(Run)),
                            by = "feature"]
    count_fractions[n_runs > 1, feature]
}


#' Remove overlapped features
#' @inheritParams .handleFractions
#' @param summary_function summary function (mean, sum, max) that will be used
#' to pick one feature from multiple overlapping features
#' @param overlapped_features features that overlap.
#' @return `data.table`
#' @keywords internal
.filterOverlapped = function(input, summary_function, overlapped_features) {
    Intensity = id = agg_intensity = max_intensity = feature =  NULL
    
    overlapped = input[feature %in% overlapped_features]
    overlapped[, agg_intensity := summary_function(Intensity, na.rm = TRUE),
               by = c("feature", "id")]
    overlapped[, max_intensity := max(agg_intensity),
               by = "feature"]
    overlapped = overlapped[agg_intensity != max_intensity]
    input[!(id %in% overlapped$id), 
          !(colnames(input) %in% c("agg_intensity", "max_intensity")), 
          with = FALSE]
}

.countCommonFeatures = function(features_1, features_2) {
    data.table::uniqueN(intersect(as.character(features_1), as.character(features_2)))
}


.checkMultiRun = function(input) {
    Run = NULL
    
    if (is.element("Fraction", colnames(input))) {
        return(list(has_fractions = TRUE, is_risky = FALSE))
    } else {
        count_techreps = input[, list(n_techreps = data.table::uniqueN(Run)),
                               by = c("Condition", "BioReplicate")]
        if (!any(count_techreps$n_techreps > 1)) {
            has_fractions = FALSE
            is_risky = FALSE
        } else {
            info = unique(input[, list(Condition, BioReplicate, Run)])
            info[, condition := paste(Condition, BioReplicate, sep = "_")]
            single_sample = unique(info[condition == unique(condition)[1],
                                        list(Condition, BioReplicate)])
            single_sample_data = input[!is.na(Intensity) & 
                                           Condition == single_sample$Condition &
                                           BioReplicate == single_sample$BioReplicate]
            single_run_features = unique(single_sample_data[Run == unique(Run)[1], 
                                                            as.character(feature)])
            common = single_sample_data[, list(n_common = .countCommonFeatures(feature, single_run_features)), 
                                        by = "Run"]
            common$fraction = common$n_common / max(common$n_common)
            overlap = common$fraction[-1]
            
            if (all(overlap > 0.5)) {
                has_fractions = FALSE
                is_risky = FALSE
            } else if (all(overlap < 0.5)) {
                has_fractions = TRUE
                is_risky = FALSE
            } else {
                has_fractions = FALSE
                is_risky = TRUE
            }
            
        }
    }
    list(has_fractions = has_fractions,
         is_risky = is_risky)
}


.handleFractionsLF = function(input) {
    check_multi_run = .checkMultiRun(input)
    
    if (check_multi_run$is_risky) {
        msg = paste("** MSstats suspects that there are fractionations and",
                    "potentially technical replicates too.",
                    "Please add Fraction column to the input.")
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    } else if (check_multi_run$has_fractions) {
        if (!is.element("Fraction", colnames(input))) {
            input = .addFractions(input)
            has_missing_fractions = any(is.na(input$Fraction))
            if (has_missing_fractions) {
                msg = paste("** It is hard to find the same fractionation across sample,",
                            "due to lots of overlapped features between fractionations.",
                            "Please add Fraction column in input.")
                getOption("MSstatsLog")("ERROR", msg)
                stop(msg)
            }
        }
        n_fractions = data.table::uniqueN(input$Fraction)
        msg = paste("Multiple fractionations exist:",
                    n_fractions, "fractionations per MS replicate.")
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        
        input = .removeOverlappingFeatures(input)
        input = .checkOverlappedFeatures(input)
        # TODO: check if there are still overlapped fractions here, stop if so?
    } else {
        input$Fraction = 1L
    }
    input
}

.addFractions = function(input) {
    input$Fraction <- NA
    run_info = unique(input[, list(Condition, BioReplicate, Run,
                                   CONDITION = paste(Condition, BioReplicate, sep = "_"))])
    single_condition = run_info[CONDITION == unique(CONDITION)[1], ]
    single_condition$Fraction = 1:nrow(single_condition)
    for (run_id in seq_along(single_condition$Run)) {
        features_in_run = as.character(unique(input[!is.na(Intensity) & Run == single_condition$Run[run_id], feature]))
        same_features = input[!is.na(Intensity) & feature %in% features_in_run, 
                              list(Condition, BioReplicate, Run, Intensity)]
        count_features = data.table::dcast(same_features, Run ~ Condition + BioReplicate,
                                           fun.aggregate = length, value.var = "Intensity")
        same_fraction = apply(count_features[, -1], 2, function(x) as.character(count_features[which.max(x)," Run"]))
        input[Run %in% same_fraction, "Fraction"] = single_condition[Run == single_condition$Run[run_id], Fraction]
    }
    input
}

.removeOverlappingFeatures = function(input) {
    fraction_keep = Fraction = NULL
    
    if (data.table::uniqueN(input$Fraction) > 1) {
        input[, fraction_keep := .getCorrectFraction(.SD), 
              by = "feature", 
              .SDcols = c("feature", "Fraction", "Run", "Intensity")] # by = c("LABEL", "PROTEINNAME", "feature")?
        input[, Intensity := ifelse(Fraction != fraction_keep, NA, Intensity)]
        input[, !(colnames(input) == "fraction_keep"), with = FALSE]
    }
    input
}


.getCorrectFraction = function(input) {
    Intensity = Run = NULL
    
    measurement_count = input[!is.na(Intensity) & Intensity > 0, 
                              list(n_obs = data.table::uniqueN(Run)),
                              by = "Fraction"]
    which_max_measurements = which.max(measurement_count$n_obs)
    if (length(which_max_measurements) == 1L) {
        return(unique(measurement_count$n_obs[which_max_measurements]))
    } else {
        average_abundance = input[!is.na(Intensity) & Intensity > 0, 
                                  list(mean_abundance = mean(Intensity)),
                                  by = "Fraction"]
        which_max_abundance = which.max(average_abundance$mean_abundance)
        unique(average_abundance$mean_abundance[which_max_abundance])
    }
}

.checkOverlappedFeatures = function(input) {
    n_fractions = Fraction = NULL 
    
    count_fractions = input[, 
                            list(n_fractions = data.table::uniqueN(Fraction)),
                            by = "feature"]
    count_fractions = count_fractions[n_fractions > 1, ]
    if (nrow(count_fractions) > 1) {
        overlapped_features = unique(as.character(count_fractions$feature))
        msg = paste("Those features are measured across all fractionations.",
                    "Please keep only one intensity of listed features", 
                    "among fractionations from one sample.\n",
                    paste(overlapped_features, sep = ", ", collapse = ", "))
        getOption("MSstatsLog")("ERROR", msg)
        stop(msg)
    } else {
        input    
    }
}
