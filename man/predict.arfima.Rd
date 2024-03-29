% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.ARFIMA.R
\name{predict.arfima}
\alias{predict.arfima}
\title{Predicts from a fitted object.}
\usage{
\method{predict}{arfima}(
  object,
  n.ahead = 1,
  prop.use = "default",
  newxreg = NULL,
  predint = 0.95,
  exact = c("default", T, F),
  setmuhat0 = FALSE,
  cpus = 1,
  trend = NULL,
  n.use = NULL,
  xreg = NULL,
  ...
)
}
\arguments{
\item{object}{A fitted \code{arfima} object}

\item{n.ahead}{The number of steps ahead to predict}

\item{prop.use}{The proportion (between 0 and 1) or percentage (between
>1 and 100) of data points to use for prediction.  Defaults to the string
"default", which sets the number of data points \code{n.use} to the minimum
of the series length and 1000.  Overriden by \code{n.use}.}

\item{newxreg}{If a regression fit, the new regressors}

\item{predint}{The percentile to use for prediction intervals assuming normal
deviations.}

\item{exact}{Controls whether exact (based on the theoretical autocovariance
matrix) prediction variances are calculated (which is recommended), as well
as whether the exact prediction formula is used when the process is
differenced (which can take a fair amount of time if the length of the series
used to predict is large).  Defaults to the string "default", which is
\code{TRUE} for the first and \code{FALSE} for the second.  A Boolean value
(\code{TRUE} or \code{FALSE}) will set both to this value.}

\item{setmuhat0}{Experimental. Sets muhat equal to zero}

\item{cpus}{The number of CPUs to use for prediction. Currently not
implemented}

\item{trend}{An optional vector the length of \code{n.ahead} or longer to
add to the predictions}

\item{n.use}{Directly set the number mentioned in \code{prop.use}.}

\item{xreg}{Alias for newxreg}

\item{\dots}{Optional arguments. Currently not used}
}
\value{
A list of lists, ceiling(prop.use * n)one for each mode with relavent details about the
prediction
}
\description{
Performs prediction of a fitted \code{arfima} object. Includes prediction
for each mode and exact and limiting prediction error standard deviations.
\strong{NOTE:  the standard errors in beta are currently not taken into
account in the prediction intervals shown.  This will be updated as soon
as possible.}
}
\examples{

\donttest{
set.seed(82365)
sim <- arfima.sim(1000, model = list(dfrac = 0.4, theta=0.9, dint = 1))
fit <- arfima(sim, order = c(0, 1, 1), back=TRUE)
fit
pred <- predict(fit, n.ahead = 5)
pred
plot(pred, numback=50)
#Predictions aren't really different due to the
#series.  Let's see what happens when we regress!

set.seed(23524)
#Forecast 5 ahead as before
#Note that we need to integrate the regressors, since time series regression
#usually assumes that regressors are of the same order as the series.
n.fore <- 5
X <- matrix(rnorm(3000+3*n.fore), ncol = 3)
X <- apply(X, 2, cumsum)
Xnew <- X[1001:1005,]
X <- X[1:1000,]
beta <- matrix(c(2, -.4, 6), ncol = 1)
simX <- sim + as.vector(X\%*\%beta)
fitX <- arfima(simX, order = c(0, 1, 1), xreg = X, back=TRUE)
fitX
#Let's compare predictions.
predX <- predict(fitX, n.ahead = n.fore, xreg = Xnew)
predX
plot(predX, numback = 50)
#With the mode we know is really there, it looks better.
fitX <- removeMode(fitX, 2)
predXnew <- predict(fitX, n.ahead = n.fore, xreg = Xnew)
predXnew
plot(predXnew, numback=50)
}

}
\references{
Veenstra, J.Q. Persistence and Antipersistence:  Theory and
Software (PhD Thesis)
}
\seealso{
\code{\link{arfima}}, \code{\link{plot.predarfima}},
\code{\link{print.predarfima}}
}
\author{
JQ (Justin) Veenstra
}
\keyword{ts}
