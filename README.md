# MSstats development version

<!-- badges: start -->
[![BioC status](http://www.bioconductor.org/shields/build/release/bioc/MSstats.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MSstats)
<!-- badges: end -->

[MSstats](https://msstats.org) is an R-based/Bioconductor package for statistical relative quantification of peptides and proteins in mass spectrometry-based proteomic experiments. 
It is applicable to multiple types of sample preparation, including label-free workflows, workflows that use stable isotope labeled reference proteins and peptides, and work-flows that use fractionation. 
It is applicable to targeted Selected Reactin Monitoring(SRM), Data-Dependent Acquisiton(DDA or shotgun), and Data-Independent Acquisition(DIA or SWATH-MS). 

## Installation 

A stable version of the package can be downloaded from Biocondutor:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MSstats")
```

or from the official Github repository:

```
source("https://install-github.me/MeenaChoi/MSstats")
```

This development version is only available on Github:

```
source("https://install-github.me/Vitek-Lab/MSstats-dev")
```

## Documentation and examples

An extensive documentation can be found on the [official website of the package](https://msstats.org).
A vignette explaining the functionalities of the package and examples can be accessed from within R:

```
vignette(package = "MSstats")
example(package = "MSstats", topic = "dataProcess") # or any other function from the package
```

## Contributing

We welcome contributions from the community. For details on how to contribute to the
development of MSsstas, please refer to the [CONTRIBUTING](https://github.com/Vitek-Lab/MSstats-dev/.github/CONTRIBUTING.md) file.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
