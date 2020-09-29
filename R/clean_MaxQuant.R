#' Clean raw output from MaxQuant
#' @param msstats_object object that inherits from MSstatsInputFiles class.
#' @param protein_id_col character, name of a column with names of proteins.
#' @param remove_by_site logical, if TRUE, proteins only identified by site will
#' be removed.
#' @param channel_columns character, regular expression that identifies channel columns
#' in TMT data.
#' @return data.table
#' @keywords internal
.cleanRawMaxQuant = function(msstats_object, protein_id_col, 
                             remove_by_site = FALSE,
                             channel_columns = "Reporterintensitycorrected"
) {
    ProteinIDs = id = PSM = PeptideSequence = NULL
    ProteingroupIDs = PrecursorCharge = Intensity = NULL
    
    mq_input = getInputFile(msstats_object, "evidence")
    mq_pg = getInputFile(msstats_object, "protein_groups")
    
    filter_cols = c("Contaminant", "Potentialcontaminant", "Reverse")
    msg = paste("** + Contaminant, + Reverse, + Potential.contaminant",
                "proteins are removed.")
    if (remove_by_site) {
        filter_cols = c(filter_cols, "Onlyidentifiedbysite")
        msg = paste("** + Contaminant, + Reverse, + Potential.contaminant,", 
                    "+ Only.identified.by.site proteins are removed.")
    }
    
    mq_input = .filterManyColumns(mq_input, filter_cols, "+")
    mq_pg = .filterManyColumns(mq_pg, filter_cols, "+")
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
    
    mq_input[, ProteingroupIDs := as.integer(as.character(ProteingroupIDs))]
    
    if (getDataType(msstats_object) == "MSstats") {
        mq_input = mq_input[ProteingroupIDs %in% unique(mq_pg[["id"]]), ]
    }
    
    mq_input = merge(mq_input, 
                     unique(mq_pg[, list(uniquefromProteinGroups = ProteinIDs,
                                         ProteingroupIDs = id)]),
                     by = "ProteingroupIDs", sort = FALSE)
    
    if (getDataType(msstats_object) == "MSstatsTMT") {
        protein_id = .findAvailable(c("Proteins", "Leadingproteins", 
                                      "Leadingrazorprotein", "Genenames"),
                                    colnames(mq_input), "Proteins")
        channels = .getChannelColumns(colnames(mq_input), channel_columns)
    } else {
        protein_id = ifelse(protein_id_col == "Proteins", 
                            "uniquefromProteinGroups",
                            "Leadingrazorprotein")
        channels = character(0)
    }
    data.table::setnames(
        mq_input, 
        c(protein_id, "Modifiedsequence", "Charge", "Rawfile"), 
        c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run"),
        skip_absent = TRUE)
    mq_input[["PeptideSequence"]] = gsub("_", "", mq_input[["PeptideSequence"]])
    mq_cols = c("ProteinName", "PeptideSequence", "Modifications", 
                "PrecursorCharge", "Run", "Intensity", 
                "Fraction", "TechReplicate", "Run", "BioReplicate",
                "PSM", "Score")
    mq_cols = intersect(c(mq_cols, channels),
                        colnames(mq_input))
    mq_input = mq_input[, mq_cols, with = FALSE]
    mq_input = unique(mq_input)
    
    if (getDataType(msstats_object) == "MSstatsTMT") {
        mq_input[, PSM := paste(PeptideSequence, PrecursorCharge, 
                                1:nrow(mq_input), sep = "_")]
        mq_input = melt(mq_input, measure.vars = channels,
                        id.vars = c("ProteinName", "PeptideSequence", 
                                    "PrecursorCharge", "PSM", "Run", "Score"),
                        variable.name = "Channel", value.name = "Intensity")
        mq_input$Channel = gsub(channel_columns, "channel", mq_input$Channel)
        mq_input$Channel = .standardizeColnames(mq_input$Channel)
        suppressWarnings({
            mq_input$Intensity = ifelse(mq_input$Intensity == 0, NA,
                                        mq_input$Intensity)
        })
        mq_input = .filterFewMeasurements(mq_input, 0, FALSE, 
                                          c("PeptideSequence", 
                                            "PrecursorCharge", "Run"))
    }
    
    mq_input = mq_input[!is.na(Intensity), ]
    .logSuccess("MaxQuant", "clean")
    mq_input
}
