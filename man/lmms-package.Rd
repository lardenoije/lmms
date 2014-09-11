\name{lmms-package}
\Rdversion{1.1}
\alias{lmms-package}
\alias{lmms}
\docType{package}
\title{
Data-driven mixed effect model splines fit and differential expression analysis
}
\description{
The package provides quality control and filtering methods as well as linear mixed effect model splines techniques for modelling and differential expression analysis to model and mine highly dimensional data sets: investNoise to identify noisy profiles and filterNoise to remove them; lmmSpline to model heterogeneous time course expression profiles; lmmsDE to perform differential expression analysis of time course data to identify differential expression over time, between   groups or time and group interaction.
}
\details{
\tabular{ll}{
Package: \tab lmms\cr
Type: \tab Package\cr
Version: \tab 1.14\cr
Date: \tab 2014-05-08\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
Functions for quality control and filtering: \code{\link[lmms]{investNoise}}, \code{\link[lmms]{filterNoise}},\code{\link[lmms]{summary.noise}},\code{\link[lmms]{plot.noise}} \cr
Functions for data modeling: \code{\link[lmms]{lmmSpline}}, \code{\link[lmms]{lmmsDE}},\code{\link[lmms]{deriv.lmmspline}},\code{\link[lmms]{predict.lmmspline}} \cr
Functions for summarization: \code{\link[lmms]{summary.lmmspline}}, \code{\link[lmms]{summary.lmmsde}} \cr
Functions for plots: \code{\link[lmms]{plot.lmmspline}}, \code{\link[lmms]{plot.lmmsde}} \cr

}
\author{
Jasmin Straube with contributions from Kim-Anh Le Cao, Emma Huang and Dominique Gorse

Maintainer: Jasmin Straube <j.straube@qfab.org>
}
\keyword{package}
