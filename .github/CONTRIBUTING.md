# Contributing to MSstats

Thank you for considering contributing to the MSstas package.
This document aims to help you take the right steps to file a bug report (report an issue), 
request a new feature, fix a typo or propose code changes.


## Code of conduct

As contributors and maintainers of MSstats project, we pledge to follow the [Carpentry Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html).

Instances of abusive, harassing or other unacceptable behavior may be reported by following our [reporting guidelines](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html#reporting-guidelines).


## Typos

Typos and grammatical errors in documentation can be fixed directly via Github. 
All changes should be made in the `.R` or `.Rmd` files rather than `.Rd` files in the `/man` directory.


## Bug reports

Bug reports can be filed in two ways:

- via the official [Google group](https://groups.google.com/forum/#!forum/msstats),

- or via [GitHub issues](https://github.com/MeenaChoi/MSstats/issues) (the preferred way). 

Before filing a bug report please make sure that:

- the problem was not solved before 
(by searching older issues and possibly other sources such as the Google group)

- you can provide a [minimal reproducible example](https://stackoverflow.com/questions/5963269/how-to-make-a-great-r-reproducible-example) so that the developers can reproduce and understand the problem. 
The [reprex](https://reprex.tidyverse.org/) package might be helpful.

To help with describing the bug, we have created an issue template which can be used after starting a new issue.


## Feature requests

GitHub issues are meant primarily for bug requests.
We suggest to make feature requests via the official [Google group](https://groups.google.com/forum/#!forum/msstats) of `MSstats`.


## Code contributions

We welcome code contributions to the `MSstats`. 
Every contribution should be made via a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) on GitHub and obey the following rules:

* if the pull request is meant to fix a bug, please start an issue earlier to describe it,
* a separate Git branch should be created for each pull request,
* the pull request must pass checks performed by the continuous integration system,
* the pull request should include relevant tests. We use [testthat](https://cran.r-project.org/package=testthat) for writing tets,
* new code should be properly documented. We use [roxygen2](https://cran.r-project.org/package=roxygen2) and
[R Markdown](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html) for documentation.  
* new code should follow the tidyverse [style guide](https://style.tidyverse.org) with the following exceptions:
    - names of functions and classes should be `camelCase`,
    - other names should use the underscore `_` and lowercase, except already existing names,
    which should not be changed to ensure backward compatibility.
* when starting a pull request, please add two reviewers:
    - [Meena Choi](https://github.com/MeenaChoi) - the main developer and maintainer of the package,
    - another one of the current maintainers (a list can be found in the DESCRIPTION file),
* before the pull request is merged, please update the NEWS.md file and version of the package, following the conventions of [semantic versioning](https://semver.org/).
