dataset = data.table::data.table(
    Run = rep(1:5, each = 2)
)
dataset_tmt = data.table::data.table(
    Run = rep(1:5, each = 2),
    Channel = rep(c("a", "b"), times = 5)
)
annotation_1 = data.table::data.table(
    Run = rep(1:5, each = 2),
    TechRep = 1:10,
    Condition = rep(c("A", "B"), times = 5)
)
annotation_2 = data.table::data.table(
    Rawfile = rep(1:4, each = 2),
    TechRep = 1:8,
    Condition = rep(c("A", "B"), times = 4)
)
annotation_3 = data.table::data.table(
    Rawfile = rep(1:5, each = 2),
    TechRep = 1:10,
    Condition = rep(c("A", "B"), times = 5),
    Channel = rep(c("a", "b"), times = 5)
)
annotation_4 = data.table::data.table(
    Rawfile = c(rep(1:4, each = 2), 5),
    TechRep = 1:9,
    Condition = c(rep(c("A", "B"), times = 4), "A"),
    Channel = c(rep(c("a", "b"), times = 4), "c")
)
# Create annotation ----
## No annotation - return NULL
tinytest::expect_null(MSstatsConvert:::.makeAnnotation(dataset, NULL))
## No additional information - return NULL
tinytest::expect_null(MSstatsConvert:::.makeAnnotation(dataset, data.table::data.table(Run = 1:5)))
## Annotation provided - return annotation
tinytest::expect_identical(MSstatsConvert:::.makeAnnotation(dataset, annotation_1),
                           annotation_1)
## Column names are updated
tinytest::expect_true(
    is.element("Run",
               colnames(MSstatsConvert:::.makeAnnotation(dataset, annotation_3,
                                                         Run = "Rawfile")))
)
# Merge annotation ----
## MSstats version: no new information in annotation
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset, NULL),
    dataset
)
## MSstats: annotation provided
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset, annotation_1),
    merge(dataset, annotation_1, sort = FALSE)
)
## MSstatsTMT: all is OK
tmt_annotation = MSstatsConvert:::.makeAnnotation(dataset, annotation_3, Run = "Rawfile")
tinytest::expect_identical(
    MSstatsConvert:::.mergeAnnotation(dataset_tmt, tmt_annotation),
    merge(dataset_tmt, tmt_annotation, sort = FALSE)
)
## MSstats: missing condition
missing_condition = MSstatsConvert:::.makeAnnotation(dataset, annotation_2, Run = "Rawfile")
tinytest::expect_error(MSstatsConvert:::.mergeAnnotation(dataset, missing_condition))
## MSstatsTMT: missing channel
missing_channel = MSstatsConvert:::.makeAnnotation(dataset, annotation_4, Run = "Rawfile")
tinytest::expect_error(MSstatsConvert:::.mergeAnnotation(dataset_tmt, missing_channel))
