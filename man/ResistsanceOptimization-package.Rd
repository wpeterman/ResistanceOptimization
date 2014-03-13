\name{ResistanceOptimization-package}
\alias{ResistanceOptimization-package}
\alias{ResistanceOptimization}
\docType{package}
\title{
About this package
}
\description{
Optimize continuous resistance surfaces using Circuitscape v4.0 or greater.

}
\details{
\tabular{ll}{
Package: \tab ResistanceOptimization\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-12\cr
License: \tab GPL-3\cr
}
This minimal package provides functions to prepare data and execute a number of functions to optimize continuous resistance surfaces using Circuitscape and nlm within R. You must have Circuitscape (4.0-Beta or higher) installed to run these functions. These methods are only applicable to continuous landscape features such as percent canopy cover, and cannot be implemented with categorical surfaces such as land cover. Output from these functions include: an AIC model rank table comparing all optimized surfaces, parameters that optimized each of the top models, coefficients from the fitted mixed effects models, bootstrap model selection summary, plots of the fitted response curves, diagnostic plots of model fit, and Circuitscape outputs for each of the optimized resistance surfaces.

***This code can currently only be used with Windows machines due its use of the Circuitscape .exe file. Anyone interested in adapting the code to accommodate command-line execution on other platforms is welcome to do so.

In order to use this package, you must have CIRCUITSCAPE installed.
The 4.0-Beta release (or 4.0 version, when released) are required:
Beta: https://github.com/Circuitscape/Circuitscape/releases
Official releases: https://code.google.com/p/circuitscape/downloads/list
}
\author{
Bill Peterman <bill.peterman@gmail.com>
}
\references{
Please cite: Peterman et al. Molecular Ecology
}



