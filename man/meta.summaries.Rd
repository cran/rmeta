\name{meta.summaries}
\title{Meta-analysis based on effect estimates}
\alias{meta.summaries}
\alias{print.meta.summaries}
\alias{plot.meta.summaries}
\alias{summary.meta.summaries}
\alias{print.summary.meta.summaries}
\description{
  Computes a summary estimate and confidence interval from a collection
  of treatment effect estimates and standard errors.  Allows fixed or 
  random effects, optional quality weights.
}
\usage{
meta.summaries(d, se, method=c("fixed", "random"), weights=NULL,
               logscale=FALSE, names=NULL, data=NULL,
               conf.level=0.95, subset=NULL,na.action=na.fail)
\method{summary}{meta.summaries}(object,conf.level=NULL,\dots)
\method{plot}{meta.summaries}(x,summary=TRUE,summlabel="Summary",
                    conf.level=NULL,colors=meta.colors(),
                    xlab=NULL,logscale=NULL,\dots)
}
\arguments{
  \item{d}{Effect estimates}
  \item{se}{standard errors for \code{d}}
  \item{method}{Standard errors and default weights from fixed or
    random-effects? }
  \item{weights}{Optional weights (eg quality weights)}
  \item{logscale}{Effect is on a log scale? (for plotting) }
  \item{names}{labels for the separate studies }
  \item{data}{optional data frame to find variables in}
  \item{conf.level}{level for confidence intervals }
  \item{subset}{Which studies to use}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to \code{\link{na.fail}}.}
  \item{x,object}{a \code{meta.summaries} object}
  \item{summary}{Plot the summary odds ratio?}
  \item{summlabel}{Label for the summary odds ratio}
  \item{colors}{see \code{\link{meta.colors}}}
  \item{xlab}{label for the effect estimate axis.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The summary estimate is a weighted average. If \code{weights} are
  specified they are used, otherwise the reciprocal of the estimated
  variance is used.

  The estimated variance is the square of \code{se} for a \code{fixed}
  analysis. For a random analysis a heterogeneity variance is estimated
  and added.

  The variance of a weighted average is a weighted average of the
  estimated variances using the squares of the weights. This is the
  square of  the summary standard error.

  With the default weights these are the standard fixed and random
  effects calculations.
}

\value{
  An object of class \code{meta.summaries}, which has
  \code{print},\code{plot},\code{summary} and \code{\link{funnelplot}}
  methods.
}
\author{Thomas Lumley}

\seealso{
  \code{\link{meta.DSL}},
  \code{\link{meta.MH}},
  \code{\link{funnelplot}},
  \code{\link{metaplot}}
}

\examples{
data(catheter)
b <- meta.DSL(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
              names=Name, subset=c(13,6,5,3,12,4,11,1,8,10,2))
d <- meta.summaries(b$logs, b$selogs, names=b$names,
                    method="random", logscale=TRUE)
}
\keyword{htest}
\keyword{hplot}
