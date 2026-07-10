#' Clean raw Diann files
#' @param msstats_object an object of class `MSstatsDIANNFiles`.
#' @inheritParams DIANNtoMSstatsFormat
#' @return data.table
#' @importFrom stats na.omit
#' @keywords internal
.cleanRawDIANN <- function(msstats_object, MBR = TRUE,
                           quantificationColumn = "FragmentQuantCorrected",
                           global_qvalue_cutoff = 0.01,
                           qvalue_cutoff = 0.01,
                           pg_qvalue_cutoff = 0.01,
                           calculateAnomalyScores = FALSE,
                           anomalyModelFeatures = c(),
                           labeledAminoAcids = NULL) {
    dn_input <- getInputFile(msstats_object, "input")
    dn_input <- data.table::as.data.table(dn_input)

    # Process quantification columns
    quantificationColumn <- .cleanDIANNProcessQuantificationColumns(dn_input, quantificationColumn)

    # Add missing columns
    dn_input <- .cleanDIANNAddMissingColumns(dn_input)

    has_channel <- !is.null(labeledAminoAcids) && "Channel" %in% colnames(dn_input)

    # Select required columns
    dn_input <- .cleanDIANNSelectRequiredColumns(dn_input, quantificationColumn, MBR,
                                                 calculateAnomalyScores,
                                                 anomalyModelFeatures,
                                                 has_channel = has_channel)

    # Split concatenated values
    dn_input <- .cleanDIANNSplitConcatenatedValues(dn_input, quantificationColumn)

    # Process fragment information
    dn_input <- .cleanDIANNProcessFragmentInfo(dn_input, quantificationColumn)

    # Clean and filter data
    dn_input <- .cleanDIANNCleanAndFilterData(dn_input, MBR, quantificationColumn,
                                              global_qvalue_cutoff,
                                              qvalue_cutoff,
                                              pg_qvalue_cutoff)

    # Rename columns
    dn_input <- .cleanDIANNRenameColumns(dn_input, quantificationColumn)

    # Assign IsotopeLabelType for protein turnover
    dn_input <- .assignDIANNIsotopeLabelType(dn_input, labeledAminoAcids, has_channel)

    .logSuccess("DIANN", "clean")
    dn_input
}

#' Process quantification columns for DIANN 2.0 format
#' @param dn_input data.table input
#' @param quantificationColumn quantification column name
#' @return updated quantification column name
#' @noRd
.cleanDIANNProcessQuantificationColumns <- function(dn_input, quantificationColumn) {
    if (quantificationColumn == "auto") {
        fragment_columns <- grep("^Fr[0-9]+Quantity$", names(dn_input), value = TRUE)
        if (length(fragment_columns) == 0) {
            stop("No fragment quantification columns found. Please check your input.")
        }
        dn_input[, FragmentQuantCorrected := do.call(paste, c(.SD, sep = ";")),
                 .SDcols = fragment_columns]
        quantificationColumn <- "FragmentQuantCorrected"
    }
    return(quantificationColumn)
}

#' Add missing required columns
#' @param dn_input data.table input
#' @return data.table with missing columns added
#' @noRd
.cleanDIANNAddMissingColumns <- function(dn_input) {
    if (!is.element("PrecursorMz", colnames(dn_input))) {
        dn_input[, PrecursorMz := NA]
    }
    if (!is.element('FragmentInfo', colnames(dn_input))) {
        dn_input[, FragmentInfo := NA]
    }
    return(dn_input)
}

#' Select required columns based on MBR setting
#' @param dn_input data.table input
#' @param quantificationColumn quantification column name
#' @param MBR logical indicating if match between runs was used
#' @return data.table with selected columns
#' @noRd
.cleanDIANNSelectRequiredColumns <- function(dn_input, quantificationColumn, MBR,
                                             calculateAnomalyScores,
                                             anomalyModelFeatures,
                                             has_channel = FALSE) {
    base_cols <- c('ProteinNames', 'StrippedSequence', 'ModifiedSequence',
                   'PrecursorCharge', quantificationColumn, 'QValue',
                   'PrecursorMz', 'FragmentInfo', 'Run')

    if (has_channel) {
        base_cols <- c(base_cols, 'Channel')
    }

    mbr_cols <- if (MBR) {
        c('LibQValue', 'LibPGQValue')
    } else {
        c('GlobalQValue', 'GlobalPGQValue')
    }

    qual_cols <- if (calculateAnomalyScores) {
        anomalyModelFeatures
    } else {
        c()
    }

    req_cols <- c(base_cols, mbr_cols, qual_cols)
    return(dn_input[, req_cols, with = FALSE])
}

#' Split concatenated values in quantification and fragment info columns
#' @param dn_input data.table input
#' @param quantificationColumn quantification column name
#' @return data.table with split values
#' @noRd
.cleanDIANNSplitConcatenatedValues <- function(dn_input, quantificationColumn) {
    split_cols <- c(quantificationColumn, "FragmentInfo")
    by_cols <- setdiff(colnames(dn_input), split_cols)
    
    dn_input <- dn_input[, lapply(.SD, function(x) unlist(tstrsplit(x, ";"))),
                         .SDcols = split_cols, 
                         by = by_cols]
    return(dn_input)
}

#' Process fragment information and add derived columns
#' @param dn_input data.table input
#' @param quantificationColumn quantification column name
#' @return data.table with processed fragment info
#' @noRd
.cleanDIANNProcessFragmentInfo <- function(dn_input, quantificationColumn) {
    # Generate fragment info if missing
    if (all(is.na(dn_input[["FragmentInfo"]]))) {
        dn_input[, FragmentInfo := paste0("Frag", 1:.N),
                 by = c("ProteinNames", "ModifiedSequence", "PrecursorCharge", "Run")]
    }
    
    # Convert quantification column to numeric
    dn_input[, (quantificationColumn) := lapply(.SD, as.numeric), 
             .SDcols = quantificationColumn]
    
    # Process fragment ion information
    dn_input[, FragmentIon := sub('\\^\\.\\*', '', FragmentInfo)]
    
    # Extract product charge
    if (any(grepl("/", dn_input$FragmentInfo))) {
        dn_input[, ProductCharge := .cleanDIANNExtractProductCharge(FragmentInfo), by = FragmentInfo]
    } else {
        dn_input[, ProductCharge := 1]
    }
    
    return(dn_input)
}

#' Extract product charge from fragment info
#' @param fragment_info fragment information string
#' @return numeric product charge
#' @noRd
.cleanDIANNExtractProductCharge <- function(fragment_info) {
    charge_part <- unlist(strsplit(fragment_info, split = "/"))[[1]]
    return(strtoi(sub("\\.\\*\\^", "", charge_part)))
}

#' Clean and filter data by removing unwanted fragments and NA values
#' @param dn_input data.table input
#' @inheritParams DIANNtoMSstatsFormat
#' @return cleaned data.table
#' @noRd
.cleanDIANNCleanAndFilterData <- function(dn_input, MBR, quantificationColumn,
                                         global_qvalue_cutoff,
                                         qvalue_cutoff, 
                                         pg_qvalue_cutoff) {
    # Remove NH3 and H2O loss fragments & remove rows with NA in quant column
    dn_input <- dn_input[!grepl("NH3|H2O", FragmentIon) & 
                             !is.na(get(quantificationColumn))]
    
    msg = paste0('** Filtering on Q.Value < ', global_qvalue_cutoff)
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    
    dn_input = dn_input[QValue >= global_qvalue_cutoff, (quantificationColumn) := 0]
    if (MBR) {
        msg = '** MBR was used to analyze the data. Now setting names and filtering'
        msg_1_mbr = paste0('-- LibPGQValue < ', pg_qvalue_cutoff)
        msg_2_mbr = paste0('-- LibQValue < ', qvalue_cutoff)
        dn_input = dn_input[LibPGQValue >= pg_qvalue_cutoff, (quantificationColumn) := 0]
        dn_input = dn_input[LibQValue >= qvalue_cutoff, (quantificationColumn) := 0]
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg_1_mbr)
        getOption("MSstatsMsg")("INFO", msg_1_mbr)
        getOption("MSstatsLog")("INFO", msg_2_mbr)
        getOption("MSstatsMsg")("INFO", msg_2_mbr)
        # getOption("MSstatsLog")("INFO", "\n")
    } else{
        msg = '** MBR was not used to analyze the data. Now setting names and filtering'
        msg_1 = paste0('-- Filtering on GlobalPGQValue < ', pg_qvalue_cutoff)
        msg_2 = paste0('-- Filtering on GlobalQValue < ', qvalue_cutoff)
        dn_input = dn_input[GlobalPGQValue >= pg_qvalue_cutoff, (quantificationColumn) := 0]
        dn_input = dn_input[GlobalQValue >= qvalue_cutoff, (quantificationColumn) := 0]
        getOption("MSstatsLog")("INFO", msg)
        getOption("MSstatsMsg")("INFO", msg)
        getOption("MSstatsLog")("INFO", msg_1)
        getOption("MSstatsMsg")("INFO", msg_1)
        getOption("MSstatsLog")("INFO", msg_2)
        getOption("MSstatsMsg")("INFO", msg_2)
        # getOption("MSstatsLog")("INFO", "\n")
    }
    
    return(dn_input)
}

#' Assign IsotopeLabelType for DIANN protein turnover workflows.
#'
#' Dispatches to one of two classification paths depending on \code{has_channel}:
#'
#' \strong{Channel-based path} (\code{has_channel = TRUE}): \code{Channel}
#' values are mapped directly to \code{IsotopeLabelType} (\code{"H"} →
#' \code{"H"}, \code{"L"} → \code{"L"}, anything else → \code{NA}), and the
#' \code{Channel} column is then dropped.  \code{labeledAminoAcids} acts solely
#' as the opt-in flag that enables this path; the amino acid codes are
#' \strong{not} used to validate or filter \code{ModifiedSequence}.
#'
#' \strong{ModifiedSequence-parsing path} (\code{has_channel = FALSE}):
#' \code{PeptideSequence} (the retained \code{ModifiedSequence}) is scanned
#' for isotope-labeled amino acids, which appear in parentheses immediately
#' after the labeled residue, in the form \code{(<label>-<aminoAcid>-H)} for
#' the heavy label or \code{(<label>-<aminoAcid>-L)} for the light label.
#' For example, \code{K(label-K-H)} marks a heavy-labeled lysine (\code{K}).
#' Here \code{<aminoAcid>} is one of the single-letter codes in
#' \code{labeledAminoAcids}, and \code{<label>} is the label name (e.g.
#' \code{SILAC} or \code{label}). Sequences with no such parenthetical
#' tags are assigned \code{IsotopeLabelType = NA}. Once classified, the
#' parenthetical annotation is stripped out of \code{PeptideSequence},
#' leaving the plain amino acid sequence.
#'
#' @param dn_input \code{data.table} after column renaming.
#' @param labeledAminoAcids Character vector of single-letter amino acid codes
#'   (e.g. \code{c("K")} or \code{c("K", "R")}), or \code{NULL} to skip
#'   protein-turnover mode entirely (backwards-compatible default).
#' @param has_channel Logical; \code{TRUE} when the raw input contained a
#'   \code{Channel} column that was retained through
#'   \code{.cleanDIANNSelectRequiredColumns}.
#' @return \code{dn_input} with \code{IsotopeLabelType} column added (and
#'   \code{Channel} removed when \code{has_channel} is \code{TRUE}).
#' @keywords internal
#' @noRd
.assignDIANNIsotopeLabelType <- function(dn_input, labeledAminoAcids, has_channel) {
    IsotopeLabelType = Channel = PeptideSequence = NULL

    if (is.null(labeledAminoAcids)) {
        return(dn_input)
    }

    if (has_channel) {
        dn_input[, IsotopeLabelType := data.table::fcase(
            Channel == "H", "H",
            Channel == "L", "L",
            default = NA_character_
        )]
        dn_input[, Channel := NULL]
    } else {
        aa_pattern <- paste0(labeledAminoAcids, collapse = "|")
        heavy_regex <- paste0("\\([^-]+-(?:", aa_pattern, ")-H\\)")
        light_regex <- paste0("\\([^-]+-(?:", aa_pattern, ")-L\\)")
        strip_regex <- paste0("\\([^-]+-(?:", aa_pattern, ")-[HL]\\)")

        dn_input <- .classifyIsotopeLabelType(dn_input, heavy_regex, light_regex)
        dn_input[, PeptideSequence := gsub(strip_regex, "", PeptideSequence, perl = TRUE)]
    }

    if (all(is.na(dn_input[["IsotopeLabelType"]]))) {
        warning("labeledAminoAcids was provided but no rows were classified as H or L. ",
                "Check that the input contains either a Channel column with H/L values ",
                "or ModifiedSequence entries with isotope-labeled amino acids", 
                "which appear in parentheses immediately after the labeled residue", 
                "e.g. (SILAC-K-H).")
    }

    dn_input
}


#' Rename columns to standardized names
#' @param dn_input data.table input
#' @param quantificationColumn quantification column name
#' @return data.table with renamed columns
#' @noRd
.cleanDIANNRenameColumns <- function(dn_input, quantificationColumn) {
    old_names <- c('ProteinNames', 'StrippedSequence', 'ModifiedSequence',
                   'PrecursorCharge', quantificationColumn, 'QValue', 
                   'PrecursorMz', 'FragmentIon', 'Run', 'ProductCharge')
    
    new_names <- c('ProteinName', 'PeptideSequence', 'PeptideModifiedSequence',
                   'PrecursorCharge', 'Intensity', 'DetectionQValue', 
                   'PrecursorMz', 'FragmentIon', 'Run', 'ProductCharge')
    
    data.table::setnames(dn_input, old = old_names, new = new_names, skip_absent = TRUE)
    
    # Clean up peptide sequence columns
    dn_input[, PeptideSequence := NULL]
    setnames(dn_input, "PeptideModifiedSequence", "PeptideSequence")
    
    return(dn_input)
}
