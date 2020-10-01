# MSstats development version

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/Vitek-Lab/MSstatsConvert.svg?branch=master)](https://travis-ci.org/Vitek-Lab/MSstatsConvert)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstatsConvert/branch/master/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstatsConvert?branch=master)
<!-- badges: end -->

[MSstats](https://msstats.org) is an R-based/Bioconductor package for statistical relative quantification of peptides and proteins in mass spectrometry-based proteomic experiments. 
It is applicable to multiple types of sample preparation, including label-free workflows, workflows that use stable isotope labeled reference proteins and peptides, and workflows that use fractionation. 
It is applicable to targeted Selected Reactin Monitoring(SRM), Data-Dependent Acquisiton(DDA or shotgun), and Data-Independent Acquisition(DIA or SWATH-MS). 
This package implements converter functions that are used by MSstats to import data from various signal processing tools.

## Installation 

This development version is only available on Github:

```
devtools::install_github("Vitek-Lab/MSstatsConvert", build_vignettes = TRUE)
```

## Documentation and examples

An extensive documentation can be found on the [official website of the package](https://msstats.org).


## Contributing

We welcome contributions from the community. For details on how to contribute to the
development of MSsstas, please refer to the [CONTRIBUTING](https://github.com/Vitek-Lab/MSstatsConvert/.github/CONTRIBUTING.md) file.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)
