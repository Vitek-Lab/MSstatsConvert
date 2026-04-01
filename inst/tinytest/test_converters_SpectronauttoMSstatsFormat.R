# Test SpectronauttoMSstatsFormat ---------------------------
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
output = SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
expect_equal(ncol(output), 11)
expect_equal(nrow(output), 372)
expect_true("Run" %in% colnames(output))
expect_true("ProteinName" %in% colnames(output))
expect_true("PeptideSequence" %in% colnames(output))
expect_true("PrecursorCharge" %in% colnames(output))
expect_true("Intensity" %in% colnames(output))
expect_true("FragmentIon" %in% colnames(output))
expect_true("ProductCharge" %in% colnames(output))
expect_true("IsotopeLabelType" %in% colnames(output))
expect_true("Condition" %in% colnames(output))
expect_true("BioReplicate" %in% colnames(output))
expect_true("Fraction" %in% colnames(output))


# Test SpectronauttoMSstatsFormat Missing Columns ---------------------------
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$F.ExcludedFromQuantification = NULL
expect_silent(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$F.FrgLossType = NULL
expect_silent(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$PG.ProteinGroups = NULL
expect_silent(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE)
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$PG.ProteinGroups = NULL
spectronaut_raw$PG.ProteinAccessions = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: PGProteinGroups"
)

spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$FG.Charge = NULL
expect_error(
    SpectronauttoMSstatsFormat(spectronaut_raw, use_log_file = FALSE),
    "The following columns are missing from the input data: FGCharge"
)


# Test SpectronauttoMSstatsFormat Quality Metrics ---------------------------
spectronaut_raw = system.file("tinytest/raw_data/Spectronaut/spectronaut_input.csv",
                              package = "MSstatsConvert")
spectronaut_raw = data.table::fread(spectronaut_raw)
spectronaut_raw$`FG.ShapeQualityScore (MS2)` = 1
spectronaut_raw$`FG.ShapeQualityScore (MS1)` = 1
spectronaut_raw$`EG.ApexRT` = 1
spectronaut_raw$`F.PossibleInterference` = TRUE
temporal = spectronaut_raw[, .(R.Replicate, R.Condition), by = R.FileName]
temporal = unique(temporal)
data.table::setnames(temporal, "R.FileName", "Run")
temporal$Order = seq(1:nrow(temporal))
temporal = temporal[, c("Run", "Order")]

msstats_format = SpectronauttoMSstatsFormat(
    spectronaut_raw,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
)

expect_true("Run" %in% colnames(msstats_format))
expect_true("ProteinName" %in% colnames(msstats_format))
expect_true("PeptideSequence" %in% colnames(msstats_format))
expect_true("PrecursorCharge" %in% colnames(msstats_format))
expect_true("Intensity" %in% colnames(msstats_format))
expect_true("FragmentIon" %in% colnames(msstats_format))
expect_true("ProductCharge" %in% colnames(msstats_format))
expect_true("IsotopeLabelType" %in% colnames(msstats_format))
expect_true("Condition" %in% colnames(msstats_format))
expect_true("BioReplicate" %in% colnames(msstats_format))
expect_true("Fraction" %in% colnames(msstats_format))
expect_true("AnomalyScores" %in% colnames(msstats_format))
expect_true("FGShapeQualityScore(MS2)" %in% colnames(msstats_format))
expect_true("FGShapeQualityScore(MS1)" %in% colnames(msstats_format))
expect_true("EGApexRT" %in% colnames(msstats_format))
expect_true("FGShapeQualityScore(MS2).mean_decrease" %in% colnames(msstats_format))
expect_true("FGShapeQualityScore(MS1).mean_increase" %in% colnames(msstats_format))
expect_true("EGApexRT.dispersion_increase" %in% colnames(msstats_format))


# Spectronaut quality features error handling -------------------------------

# runOrder not provided
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = NULL,
    n_trees = 100,
    max_depth = "auto"
))

# runOrder is not a data frame
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = seq(1:nrow(annotation)),
    n_trees = 100,
    max_depth = "auto"
))

# n_trees is negative
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = -3,
    max_depth = "auto"
))

# FAIL n_trees is a decimal
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100.5,
    max_depth = "auto"
))

# n_trees is a string
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = "string",
    max_depth = "auto"
))

# max_depth is a string
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "string"
))

# FAIL max_depth is a decimal
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = 5.5
))

# FAIL max_depth is a negative number
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c("FG.ShapeQualityScore (MS2)",
                             "FG.ShapeQualityScore (MS1)",
                             "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = -5
))

# anomalyModelFeatures is empty
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures = c(),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# anomalyModelFeatureTemporal contains unrecognized string
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("string",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# anomalyModelFeatureTemporal is empty
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c(),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# FAIL removeMissingFeatures is a percentage above 1
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 1.5,
    anomalyModelFeatureCount = 100,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# FAIL anomalyModelFeatureCount is not an integer
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = 100.5,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# FAIL anomalyModelFeatureCount is a string
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = "string",
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))

# anomalyModelFeatureCount is a negative number
expect_error(SpectronauttoMSstatsFormat(
    spectronaut_raw,
    annotation = annotation,
    calculateAnomalyScores = TRUE,
    anomalyModelFeatures =c("FG.ShapeQualityScore (MS2)",
                            "FG.ShapeQualityScore (MS1)",
                            "EG.ApexRT"),
    anomalyModelFeatureTemporal = c("mean_decrease",
                                    "mean_increase",
                                    "dispersion_increase"),
    removeMissingFeatures = 0.5,
    anomalyModelFeatureCount = -5,
    runOrder = temporal,
    n_trees = 100,
    max_depth = "auto"
))
