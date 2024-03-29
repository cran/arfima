% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arfima_function.R
\name{weed}
\alias{weed}
\title{Weeds out fits from a call to arfima that are too close to each other.}
\usage{
weed(
  ans,
  type = c("A", "P", "B", "N"),
  walls = FALSE,
  eps2 = 0.025,
  eps3 = 0.01,
  adapt = TRUE,
  pn = 2
)
}
\arguments{
\item{ans}{The result from a call to arfima.}

\item{type}{The space to perform the weeding in. "A" is for operating
parameters. "P" is in the PACF space. "B" performs weeding in both. "N"
performs no weeding and is only used internally.}

\item{walls}{If more than one mode is on a wall in the PACF space, all modes
but the one with the highest log-likelihood on the same wall are deleted.}

\item{eps2}{The maximum distance between modes that are close together for
the mode with the lower log-likelihood to be weeded out.  If adapt is TRUE
(default) this value changes.}

\item{eps3}{The minimum distance from a wall for a secondary mode to be
weeded out, if walls are TRUE.}

\item{adapt}{If TRUE, if dim is the dimensionality of the search, eps2 is
changed to \eqn{(1 + eps2)^{dim} - 1}.}

\item{pn}{The p in the p-norm to be used in the weeding.  p = 2 (default) is
Euclidean distance.}
}
\value{
An object of class "arfima" with modes possibly weeded out.
}
\description{
Weeds out fits from a call to arfima that are too close to each other.
}
\examples{

\donttest{
set.seed(1234)
sim <- arfima.sim(1000, model = list(theta = 0.9, dfrac = 0.4))
fit <- arfima(sim, order = c(0, 0, 1), autoweed = FALSE, back=TRUE)
fit
distance(fit)
fit1 <- weed(fit)
fit1
distance(fit1)
}
}
\seealso{
\code{\link{arfima}}, \code{\link{distance}}
}
\author{
JQ (Justin) Veenstra
}
\keyword{ts}
