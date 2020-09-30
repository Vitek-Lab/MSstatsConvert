dataset = data.table::data.table(
    Run = 1:5,
    Condition = 1,
    BioReplicate = 1
)
dataset_tmt = data.table::data.table(
    Run = as.character(rep(1:5, each = 2)),
    Channel = rep(c("a", "b"), times = 5)
)
annotation_1 = data.table::data.table(
    Run = 1:5,
    Condition = 1,
    BioReplicate = 1
)
annotation_2 = data.table::data.table(
    Rawfile = 1:4,
    Condition = 1,
    BioReplicate = 1
)
annotation_3 = data.table::data.table(
    Rawfile = rep(1:5, each = 2),
    Condition = rep(c("A", "B"), times = 5),
    Channel = rep(c("a", "b"), times = 5),
    Mixture = 1,
    TechRepMixture = 1,
    Fraction = 1,
    BioReplicate = 1
)
annotation_4 = data.table::data.table(
    Rawfile = c(rep(1:4, each = 2), 5),
    TechRep = 1:9,
    Condition = c(rep(c("A", "B"), times = 4), "A"),
    Channel = c(rep(c("a", "b"), times = 4), "c"),
    Mixture = 1,
    TechRepMixture = 1,
    Fraction = 1,
    BioReplicate = 1
)
# Create annotation ----
## No annotation - return a subset of original data
unique_runs = unique(dataset)
unique_runs$Run = as.character(unique_runs$Run)
tinytest::expect_equal(MSstatsConvert:::MSstatsMakeAnnotation(dataset, NULL),
                       unique_runs)
## No additional information - return NULL
tinytest::expect_equal(MSstatsConvert:::MSstatsMakeAnnotation(
    dataset, data.table::data.table(Run = 1:5,
                                    Condition = 1,
                                    BioReplicate = 1)),
    data.table::data.table(Run = as.character(1:5),
                           Condition = 1,
                           BioReplicate = 1))
## Annotation provided - return annotation
annotation_1$Run = as.character(annotation_1$Run)
tinytest::expect_identical(MSstatsConvert:::MSstatsMakeAnnotation(dataset, annotation_1),
                           annotation_1)
## Column names are updated
tinytest::expect_true(
    is.element("Run",
               colnames(MSstatsConvert:::MSstatsMakeAnnotation(dataset, 
                                                               annotation_2,
                                                               Run = "Rawfile")))
)
# Merge annotation ----
## MSstats version: no new information in annotation
dataset_chr = dataset
dataset_chr$Run = as.character(dataset_chr$Run)
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset, NULL),
    dataset_chr
)
## MSstats: annotation provided
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset, annotation_1),
    merge(dataset_chr, annotation_1, sort = FALSE)
)
## MSstatsTMT: all is OK
tmt_annotation = MSstatsConvert:::MSstatsMakeAnnotation(dataset_tmt, annotation_3, Run = "Rawfile")
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset_tmt, tmt_annotation),
    merge(dataset_tmt, tmt_annotation, sort = FALSE)
)
## MSstats: missing condition
missing_condition = MSstatsConvert:::MSstatsMakeAnnotation(dataset, annotation_2, Run = "Rawfile")
tinytest::expect_message(MSstatsConvert:::.mergeAnnotation(dataset, missing_condition))
## MSstatsTMT: missing channel
missing_channel = MSstatsConvert:::MSstatsMakeAnnotation(dataset_tmt, annotation_4, Run = "Rawfile")
tinytest::expect_error(MSstatsConvert:::.mergeAnnotation(dataset_tmt, missing_channel))

