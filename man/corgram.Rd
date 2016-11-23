\name{corgram}
\title{Calculate and Plot a Correlogram}
\alias{corgram}
\usage{
corgram(y, wt=1, maxlag=NULL, partial=FALSE, add=FALSE, lty=1, xlim=NULL,
	ylim=NULL, xlab=NULL, ylab=NULL, main=NULL, ...)
}
\description{
\code{corgram} calculates the values of a correlogram (autocorrelation
function or ACF) and plots it.
}
\arguments{
\item{y}{A time series vector.}
\item{maxlag}{Maximum number of lags for which the correlation is to
be calculated.}
\item{partial}{If TRUE, the partial autocorrelation function (PACF) is
plotted.}
\item{wt}{Indicator vector with zeros for values to be ignored.}
\item{add}{If TRUE, adds a new correlogram to an existing plot.}
\item{others}{Plotting parameters}
}
\value{
\code{corgram} returns a two-column matrix containing the (partial) correlogram
coordinates.
}
\author{J.K. Lindsey}
\examples{
y <- rnorm(100)
corgram(y)
corgram(y, partial=TRUE)
}
\keyword{hplot}
