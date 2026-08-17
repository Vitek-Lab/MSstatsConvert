# Comprehensive Unit Tests for Anomaly Detection Function using tinytest

# Helper function to create test data with quality metrics
run_quality_metrics = function(df, mean_increase, mean_decrease, dispersion_increase) {
    df$QualityMetric.mean_increase = mean_increase
    df$QualityMetric.mean_decrease = mean_decrease
    df$QualityMetric.dispersion_increase = dispersion_increase
    anomaly_output = MSstatsConvert:::.runAnomalyModel(df, 
                                                        n_trees=100, 
                                                        max_depth="auto", 
                                                        cores=1,
                                                        split_column="PSM",
                                                        quality_metrics=c("QualityMetric.mean_increase",
                                                                          "QualityMetric.mean_decrease",
                                                                          "QualityMetric.dispersion_increase"))
    return(anomaly_output)
}

# Create base dataframe for testing
create_base_df = function(n_rows) {
    data.table::data.table(
        PSM = rep("PSM1", n_rows)
    )
}

# Test 1: Property-Based Testing - Monotonicity with Cumulative Sums
base_df_10 = create_base_df(10)

# Baseline data with low values
baseline_scores = run_quality_metrics(
    base_df_10,
    rep(0.1, 10),  # mean_increase
    rep(0.1, 10),  # mean_decrease  
    rep(0.1, 10)   # dispersion_increase
)

# Data with progressively higher cumulative sums
high_scores = run_quality_metrics(
    base_df_10,
    c(seq(0, 0.1, length.out = 5), seq(2.0, 5.0, length.out = 5)),  # mean_increase
    c(seq(0, 0.1, length.out = 5), seq(2.0, 5.0, length.out = 5)),  # mean_decrease
    c(seq(0, 0.1, length.out = 5), seq(2.0, 5.0, length.out = 5))   # dispersion_increase
)

# The last 5 rows (with high values) should have higher mean anomaly scores
expect_true(mean(high_scores$AnomalyScores[6:10]) > mean(high_scores$AnomalyScores[1:5]),
            info = "Higher cumulative sum values should produce higher anomaly scores")

# Test 2: Extreme Value Testing - Obvious Outliers
base_df_20 = create_base_df(20)

extreme_scores = run_quality_metrics(
    base_df_20,
    c(seq(0, 0.1, length.out = 19), 10.0),  # Last value is extreme
    c(seq(0, 0.1, length.out = 19), 8.0),   # Last value is extreme
    c(seq(0, 0.1, length.out = 19), 12.0)   # Last value is extreme
)

# The extreme outlier (last row) should have the highest anomaly score
expect_true(extreme_scores$AnomalyScores[20] == max(extreme_scores$AnomalyScores),
            info = "Extreme outlier should have highest anomaly score")

# The outlier should score significantly higher than the median
expect_true(extreme_scores$AnomalyScores[20] > median(extreme_scores$AnomalyScores[1:19]) * 2,
            info = "Outlier should score significantly higher than median")

# Test 3: Consistency/Reproducibility Testing
base_df_20_orig = create_base_df(20)

# Run multiple times with original test data
scores1 = run_quality_metrics(
    base_df_20_orig,
    c(0.00, 0.00, 0.26, 0.00, 0.00, 
      0.00, 0.50, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.55, 
      0.87, 1.11, 1.42, 1.71, 2.94),
    c(0.00, 0.00, 0.00, 0.29, 0.50, 
      0.56, 0.00, 1.16, 0.99, 0.15, 
      0.00, 0.00, 1.27, 1.84, 0.29, 
      0.00, 0.00, 0.00, 0.00, 0.00),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.01, 0.85, 0.00, 0.00, 
      0.00, 0.00, 0.96, 1.07, 1.15, 
      0.89, 0.50, 0.22, 0.00, 0.91)
)

scores2 = run_quality_metrics(
    base_df_20_orig,
    c(0.00, 0.00, 0.26, 0.00, 0.00, 
      0.00, 0.50, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.55, 
      0.87, 1.11, 1.42, 1.71, 2.94),
    c(0.00, 0.00, 0.00, 0.29, 0.50, 
      0.56, 0.00, 1.16, 0.99, 0.15, 
      0.00, 0.00, 1.27, 1.84, 0.29, 
      0.00, 0.00, 0.00, 0.00, 0.00),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.01, 0.85, 0.00, 0.00, 
      0.00, 0.00, 0.96, 1.07, 1.15, 
      0.89, 0.50, 0.22, 0.00, 0.91)
)

# Correlation between runs should be high (>0.7) even if not identical
correlation = cor(scores1$AnomalyScores, scores2$AnomalyScores)
expect_true(correlation > 0.95,
            info = "Multiple runs should produce correlated results")

# Top anomalies should be largely consistent
top5_run1 = order(scores1$AnomalyScores, decreasing = TRUE)[1:5]
top5_run2 = order(scores2$AnomalyScores, decreasing = TRUE)[1:5]
overlap = length(intersect(top5_run1, top5_run2))
expect_true(overlap >= 3,
            info = "At least 3 of top 5 anomalies should overlap between runs")

# Test 4: Expected Anomaly Identification from Original Data
original_anomaly_scores = run_quality_metrics(
    base_df_20_orig,
    c(0.00, 0.00, 0.26, 0.00, 0.00, 
      0.00, 0.50, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.55, 
      0.87, 1.11, 1.42, 1.71, 2.94),
    c(0.00, 0.00, 0.00, 0.29, 0.50, 
      0.56, 0.00, 1.16, 0.99, 0.15, 
      0.00, 0.00, 1.27, 1.84, 0.29, 
      0.00, 0.00, 0.00, 0.00, 0.00),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.01, 0.85, 0.00, 0.00, 
      0.00, 0.00, 0.96, 1.07, 1.15, 
      0.89, 0.50, 0.22, 0.00, 0.91)
)

# Row 20 has highest mean_increase (2.94) - should be highly anomalous
top_anomalies = order(original_anomaly_scores$AnomalyScores, decreasing = TRUE)[1:5]
expect_true(20 %in% top_anomalies,
            info = "Row 20 with highest mean_increase should be in top anomalies")

# Create test data vectors for comparison
mean_increase_vals = c(0.00, 0.00, 0.26, 0.00, 0.00, 
                        0.00, 0.50, 0.00, 0.00, 0.00, 
                        0.00, 0.00, 0.00, 0.00, 0.55, 
                        0.87, 1.11, 1.42, 1.71, 2.94)
mean_decrease_vals = c(0.00, 0.00, 0.00, 0.29, 0.50, 
                        0.56, 0.00, 1.16, 0.99, 0.15, 
                        0.00, 0.00, 1.27, 1.84, 0.29, 
                        0.00, 0.00, 0.00, 0.00, 0.00)
dispersion_vals = c(0.00, 0.00, 0.00, 0.00, 0.00, 
                     0.00, 0.01, 0.85, 0.00, 0.00, 
                     0.00, 0.00, 0.96, 1.07, 1.15, 
                     0.89, 0.50, 0.22, 0.00, 0.91)

# Rows with all zeros should score lower than rows with high values
zero_rows = which(mean_increase_vals == 0 & mean_decrease_vals == 0 & dispersion_vals == 0)
high_value_rows = c(14, 16, 17, 18, 19, 20)  # Rows with notably high values

expect_true(mean(original_anomaly_scores$AnomalyScores[high_value_rows]) > 
                mean(original_anomaly_scores$AnomalyScores[zero_rows]),
            info = "High value rows should score higher than zero rows")

# Test 5: Edge Cases and Boundary Conditions
base_df_10_edge = create_base_df(10)

# Test with all identical values
identical_scores = run_quality_metrics(
    base_df_10_edge,
    rep(1.0, 10),  # All identical
    rep(1.0, 10),  # All identical
    rep(1.0, 10)   # All identical
)

expect_true(is.data.frame(identical_scores) && nrow(identical_scores) == 10,
            info = "Function should handle identical values gracefully")

# Test with minimum viable dataset
base_df_3 = create_base_df(3)
minimal_scores = run_quality_metrics(
    base_df_3,
    c(0.1, 0.2, 5.0),  # Clear outlier
    c(0.1, 0.2, 0.1),
    c(0.1, 0.1, 0.1)
)

# The outlier (row 3) should have highest score
expect_equal(which.max(minimal_scores$AnomalyScores), 3,
             info = "Clear outlier should have highest score in minimal dataset")

# Test 6: Scale Invariance Testing
base_df_5 = create_base_df(5)

# Original data
original_scale_scores = run_quality_metrics(
    base_df_5,
    c(0.1, 0.2, 0.3, 1.0, 2.0),
    c(0.1, 0.1, 0.5, 0.2, 0.1),
    c(0.1, 0.3, 0.1, 0.8, 0.2)
)

# Scaled data (multiply by 10)
scaled_scale_scores = run_quality_metrics(
    base_df_5,
    c(0.1, 0.2, 0.3, 1.0, 2.0) * 10,
    c(0.1, 0.1, 0.5, 0.2, 0.1) * 10,
    c(0.1, 0.3, 0.1, 0.8, 0.2) * 10
)

# Relative ranking should be preserved under scaling
original_ranking = order(original_scale_scores$AnomalyScores, decreasing = TRUE)
scaled_ranking = order(scaled_scale_scores$AnomalyScores, decreasing = TRUE)

# Rankings should be identical or very similar
expect_true(cor(original_ranking, scaled_ranking, method = "spearman") > 0.95,
            info = "Relative rankings should be preserved under scaling")

# Test 7: Regression Testing with Original Data
regression_scores = run_quality_metrics(
    base_df_20_orig,
    c(0.00, 0.00, 0.26, 0.00, 0.00, 
      0.00, 0.50, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.55, 
      0.87, 1.11, 1.42, 1.71, 2.94),
    c(0.00, 0.00, 0.00, 0.29, 0.50, 
      0.56, 0.00, 1.16, 0.99, 0.15, 
      0.00, 0.00, 1.27, 1.84, 0.29, 
      0.00, 0.00, 0.00, 0.00, 0.00),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.01, 0.85, 0.00, 0.00, 
      0.00, 0.00, 0.96, 1.07, 1.15, 
      0.89, 0.50, 0.22, 0.00, 0.91)
)

# Basic sanity checks
expect_equal(nrow(regression_scores), 20,
             info = "Output should have same number of rows as input")
expect_true(all(regression_scores$AnomalyScores >= 0),
            info = "All anomaly scores should be non-negative")
expect_true(all(is.finite(regression_scores$AnomalyScores)),
            info = "All anomaly scores should be finite")

# Row 20 (highest mean_increase = 2.94) should be highly anomalous
expect_true(regression_scores$AnomalyScores[20] > quantile(regression_scores$AnomalyScores, 0.8),
            info = "Row with highest mean_increase should be highly anomalous")

# Store results for regression testing (compare against future runs)
cat("Top 5 anomalous rows:", order(regression_scores$AnomalyScores, decreasing = TRUE)[1:5], "\n")
cat("Score range:", range(regression_scores$AnomalyScores), "\n")

# Test 8: Zero/Null Data Handling
base_df_10_zero = create_base_df(10)

# All zeros data
zero_scores = run_quality_metrics(
    base_df_10_zero,
    rep(0, 10),  # All zeros
    rep(0, 10),  # All zeros
    rep(0, 10)   # All zeros
)

expect_true(is.data.frame(zero_scores) && nrow(zero_scores) == 10,
            info = "Function should handle all-zero data gracefully")

# All scores should be similar/low for identical data
expect_true(sd(zero_scores$AnomalyScores) < 0.1,
            info = "All-zero data should have low variance in scores")

# Test 9: Relative Ranking Preservation
base_df_6_rank = create_base_df(6)

# Create data with obvious ranking: Row 6 > Row 5 > Row 4 > Rows 1,2,3
ranking_scores = run_quality_metrics(
    base_df_6_rank,
    c(0.1, 0.11, 0.12, 1.0, 2.0, 5.0),
    c(0.1, 0.11, 0.12, 1.0, 2.0, 5.0),
    c(0.1, 0.11, 0.12, 1.0, 2.0, 5.0)
)

# Row 5 should have highest score, Row 4 second highest, etc.
expect_true(ranking_scores$AnomalyScores[6] > ranking_scores$AnomalyScores[5],
            info = "Row 6 should score higher than Row 5")
expect_true(ranking_scores$AnomalyScores[5] > ranking_scores$AnomalyScores[4],
            info = "Row 5 should score higher than Row 4")
expect_true(ranking_scores$AnomalyScores[4] > max(ranking_scores$AnomalyScores[1:3]),
            info = "Row 4 should score higher than Rows 1-3")

# Test 10: Original Quality Metrics Calculation Test (from the beginning of the file)
# Test add_increase, add_decrease, add_dispersion
quality_vector = c(
    1.31, 0.05, 0.76, -0.79, -0.71, 
    -0.56, 1.00, -1.66, -0.33,  0.34, 
    0.00, -0.40, -1.77, -1.07, 1.05,
    0.82, 0.74, 0.81, 0.79, 1.73
)

mean_increase = MSstatsConvert:::.add_mean_increase(quality_vector)
mean_decrease = MSstatsConvert:::.add_mean_decrease(quality_vector)
dispersion_increase = MSstatsConvert:::.add_dispersion_increase(quality_vector)

expect_equal(
    mean_increase,
    c(0.00, 0.00, 0.26, 0.00, 0.00, 
      0.00, 0.50, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.55, 
      0.87, 1.11, 1.42, 1.71, 2.94),
    info = "Mean increase calculation should match expected values"
)

expect_equal(
    mean_decrease,
    c(0.00, 0.00, 0.00, 0.29, 0.50, 
      0.56, 0.00, 1.16, 0.99, 0.15, 
      0.00, 0.00, 1.27, 1.84, 0.29, 
      0.00, 0.00, 0.00, 0.00, 0.00),
    info = "Mean decrease calculation should match expected values"
)

expect_equal(
    dispersion_increase,
    c(0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.01, 0.85, 0.00, 0.00, 
      0.00, 0.00, 0.96, 1.07, 1.15, 
      0.89, 0.50, 0.22, 0.00, 0.91),
    tolerance = 1e-2,
    info = "Dispersion increase calculation should match expected values"
)

cat("All anomaly detection tests completed successfully!\n")

# Test n_feat and missing_run_count parameters in .prepareSpectronautAnomalyInput
spectronaut_raw = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, NA, 2200, NA, NA,
                  1300, 1800, NA, 2800, NA,
                  100, 200, 300, 400, 500),
    QualityMetric1 = rep(0.5, 25)
)

temporal = spectronaut_raw[, .(Condition = unique(Condition), 
                                    Run = unique(Run)), 
                                by = BioReplicate]
temporal$Order = seq(1:nrow(temporal))
temporal = temporal[, c("Run", "Order")]

missing_run_excluded = MSstatsConvert:::.prepareSpectronautAnomalyInput(
    spectronaut_raw, 
    quality_metrics = c("QualityMetric1"),
    run_order = temporal,
    n_feat = 2,
    missing_run_count = 0.5)
expect_false("AFPLAEWQPSDVDQR" %in% missing_run_excluded$PeptideSequence)
expect_true("LowAbundancePeptide" %in% missing_run_excluded$PeptideSequence)

low_abundance_excluded = MSstatsConvert:::.prepareSpectronautAnomalyInput(
    spectronaut_raw, 
    quality_metrics = c("QualityMetric1"),
    run_order = temporal,
    n_feat = 2,
    missing_run_count = 0.95)
expect_true("AFPLAEWQPSDVDQR" %in% low_abundance_excluded$PeptideSequence)
expect_false("LowAbundancePeptide" %in% low_abundance_excluded$PeptideSequence)


# Test 11: Testing duplicity of quality metrics, applicable considering
# multiple fragments share the same precursor level metrics

# Data with progressively higher cumulative sums
duplicate_metrics = run_quality_metrics(
    base_df_10,
    c(rep(0.1, 5), seq(2.0, 4.0, length.out = 5)),  # mean_increase
    c(rep(0.1, 5), seq(2.0, 4.0, length.out = 5)),  # mean_decrease
    c(rep(0.1, 5), seq(2.0, 4.0, length.out = 5))   # dispersion_increase
)

# The last 5 rows (with high values) should have lower mean anomaly scores
# Since they are all clumped between 2 and 4, whereas 0.1 is by itself
expect_true(mean(duplicate_metrics$AnomalyScores[6:10]) < mean(duplicate_metrics$AnomalyScores[1:5]),
            info = "Rows 6-10 (values clumped 2-4) should have lower
            anomaly scores than rows 1-5 (isolated value of 0.1)")

nan_first_row_df = create_base_df(5)
nan_first_row_df$QualityMetric.mean_increase = c(NA, 0.2, 0.4, 0.6, 0.8)

nan_first_row_result = tryCatch({
    MSstatsConvert:::.runAnomalyModel(
        nan_first_row_df,
        n_trees = 100,
        max_depth = "auto",
        cores = 1,
        split_column = "PSM",
        quality_metrics = c("QualityMetric.mean_increase"))
}, error = function(e) e)

expect_false(inherits(nan_first_row_result, "error"),
    info = paste(
        "Anomaly model should not crash/error when a quality metric has a",
        "leading NA/NaN value within a PSM group.",
        if (inherits(nan_first_row_result, "error"))
            paste("Got error:", conditionMessage(nan_first_row_result)) else ""))

if (!inherits(nan_first_row_result, "error")) {
    expect_true(all(is.finite(nan_first_row_result$AnomalyScores)),
        info = "Anomaly scores should be finite even when a quality metric has a leading missing value")
}
