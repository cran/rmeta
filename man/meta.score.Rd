\name{meta.score}
\alias{meta.score}
\alias{print.meta.score}
\alias{summary.meta.score}
\alias{print.summary.meta.score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fixed effects (likelihood score) meta-analysis}
\description{
  Computes individual relative risks, a summary relative risk, and
  confidence intervals.  The confidence intervals are based on
  inverting the score test (similar to the so-called "Cornfield" method
  for odd ratios).
}
\usage{
meta.score(ntrt, nctrl, ptrt, pctrl, conf.level = 0.95,
names = NULL, data = NULL, subset = NULL, na.action = na.fail,
 statistic = "RR", control = 1.5)
\method{summary}{meta.score}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntrt}{Number of subjects in treated/exposed group}
  \item{nctrl}{Number of subjects in control group }
  \item{ptrt}{Number of events in treated/exposed group}
  \item{pctrl}{Number of events in control group}
  \item{names}{names or labels for studies}
  \item{data}{data frame to interpret variables}
  \item{subset}{subset of studies to include}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to \code{\link{na.fail}}.}
  \item{statistic}{"RR" for relative risk is currently the only option}
  \item{control}{used to set search limits in finding interval
    endpoints}
  \item{conf.level}{Confidence level for intervals}
  \item{object}{\code{meta.score} object}
  \item{...}{For future expansion}
}


\value{
  An object of class \code{meta.MH} with \code{print},
  % \code{plot}
  % \code{funnelplot}
  and \code{summary} methods.
  
}
\references{Gart JJ, Nam J. (1988) "Approximate Interval Estimation of
  the Ratio of Binomial Parameters: A review and corrections for
  skewness". Biometrics 44(2) 323--338}

\note{The correction for skewness described in the reference is not yet implemented.}


\seealso{\code{\link{meta.MH}}, \code{\link{metaplot}}  }
\examples{
## From Gart and Nam
mice<-data.frame(n1=c(16, 16, 18, 15),
    n0 = c(79, 87, 90, 82), 
    x1 = c(4, 2, 4, 1),
    x0 = c(5, 3, 10, 3)
)

 summary(meta.score(ntrt=n1,nctrl=n0,ptrt=x1,pctrl=x0,data=mice))

}
\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line

