# MSstatsConvert 1.22.0

## Protein Turnover Support

* **DIA-NN**: `DIANNtoMSstatsFormat` now accepts a `labeledAminoAcids` parameter to enable protein turnover workflows. Provide the single-letter amino acid codes that carry the SILAC label (e.g. `c("K")` or `c("K", "R")`). The converter automatically distinguishes heavy and light peptide forms — either from a `Channel` column in DIA-NN 2.x exports, or by parsing SILAC modification tags (e.g. `(SILAC-K-H)`) in DIA-NN 1.x exports — and populates the `IsotopeLabelType` column required by MSstats downstream.

* **Spectronaut**: `SpectronauttoMSstatsFormat` now accepts a `heavyLabels` parameter to classify peptides in protein turnover experiments. Supply the label names as they appear in square brackets in the peptide sequence (e.g. `c("Lys6")` or `c("Lys6", "Arg10")`). Peptides are automatically assigned as heavy (`H`), light (`L`), or unlabeled. Any novel label name reported by Spectronaut (e.g. `"Leu6"`, `"Phe10"`) is supported.

* **Spectronaut**: A new `peptideSequenceColumn` parameter lets you specify which column in your Spectronaut export contains the peptide sequence. This is especially useful for protein turnover reports, which may use a different column layout than a standard Spectronaut output.

* **Spectronaut**: Protein turnover reports that lack fragment-level columns (e.g. `FFrgLossType`, `FFrgIon`) are now handled automatically — previously these required manual pre-processing before import.

* **Fraction selection**: The algorithm that selects the best fraction for each peptide feature now correctly accounts for heavy and light isotope channels in protein turnover data, ensuring that light-channel coverage is prioritized when choosing fractions. When multiple fractions have identical coverage, the fraction with the highest mean intensity is chosen.

## DIA-NN Quality Control: Anomaly Detection

* `DIANNtoMSstatsFormat` now supports automated data quality scoring via an isolation-forest anomaly detection model (`calculateAnomalyScores = TRUE`). When enabled, each precursor-level feature receives an anomaly score that can be used downstream in MSstats to down-weight measurements that look like outliers. Users supply the quality metric columns to use as model features (`anomalyModelFeatures`) and, optionally, a run order table (`runOrder`) to engineer temporal trend features that detect gradual instrument drift across an experiment.

## Spectronaut: MS1 Quantification

* `SpectronauttoMSstatsFormat` now accepts `"FG.MS1Quantity"` (and other raw Spectronaut column names) as the `intensity` argument, in addition to the existing `"PeakArea"` and `"NormalizedPeakArea"` options. 

# MSstatsConvert 1.16.0

* Improved memory management for efficiency.
* Added support for ProteinProspector data conversion to MSstatsTMT format.
* Updated NA handling for Skyline and Spectronaut converters
* Added MetaMorpheus converter functionality.

# MSstatsConvert 1.9.1

* Added DIA-NN converter.

# MSstatsConvert 1.5.1

* Fixed bugs related to data types.
* Fixed bug that resulted in missing TechReplicate column in output.
* Added Philosopher converter.
* Added methods for MSstatsValidated objects.

# MSstatsConvert 1.2.1

* Fixed a bug related to SRM inputs.

# MSstatsConvert 1.1.2

* Updated Spectronaut converter to allow annotation in input file.


# MSstatsConvert 1.1.1

* Minor bug fixes
* Updated MSstatsLogsSettings function to allow differents settings for different packages

# MSstatsConvert 0.99.0

* Added logging utitilies: `MSstatsLogsSettings`, `MSstatsSaveSesionInfo`.
* Added functions for importing and cleaning data: `MSstatsImport`, `MSstatsClean`.
* Added functions for preprocessing: `MSstatsMakeAnnotation`, `MSstatsPreprocess`,
`MSstatsBalancedDesign`.
* Added vignette.

# MSstatsConvert 0.0.1

* Added a `NEWS.md` file to track changes to the package.
