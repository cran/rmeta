\name{meta.DSL}
\title{Random effects (DerSimonian-Laird) meta-analysis}
\usage{
meta.DSL(ntrt, nctrl, ptrt, pctrl, conf.level=0.95,
         names=NULL, data=NULL, subset=NULL, na.action=na.fail,statistic="OR")
\method{summary}{meta.DSL}(object, conf.level=NULL, \dots)
\method{plot}{meta.DSL}(x, summary=TRUE, summlabel="Summary",
             conf.level=NULL, colors=meta.colors(), xlab=NULL,\dots)
}
\alias{meta.DSL}
\alias{summary.meta.DSL}
\alias{print.meta.DSL}
\alias{print.summary.meta.DSL}
\alias{plot.meta.DSL}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ntrt}{Number of subjects in treated/exposed group}
 \item{nctrl}{Number of subjects in control group }
 \item{ptrt}{Number of events in treated/exposed group}
 \item{pctrl}{Number of events in control group}
 \item{conf.level}{Coverage for confidence intervals}
 \item{names}{names or labels for studies}
 \item{data}{data frame to interpret variables}
 \item{subset}{subset of studies to include}
 \item{na.action}{a function which indicates what should happen when
   the data contain \code{NA}s.  Defaults to \code{\link{na.fail}}.}
 \item{statistic}{"OR" for odds ratio, "RR" for relative risk}
 \item{x,object}{a \code{meta.DSL} object}
 \item{summary}{Plot the summary odds ratio?}
 \item{summlabel}{Label for the summary odds ratio}
 \item{colors}{see \code{\link{meta.colors}}}
 \item{xlab}{x-axis label, default is based on \code{statistic}}
 \item{\dots}{further arguments to be passed from or to methods.}
}
\description{
Computes the individual odds ratios or relative risks, the summary, the random
effects variance, and Woolf's test for heterogeneity. The \code{print}
method gives the summary and test for heterogeneity; the \code{summary}
method also gives all the individual odds ratios and confidence
intervals. Studies with zero or infinite odds ratio are omitted, as
their variance cannot be calculated sensibly.

The \code{plot} method draws a standard meta-analysis plot. The
confidence interval for each study is given by a horizontal line, and
the point estimate is given by a square whose height is inversely
proportional to the standard error of the estimate.  The summary odds
ratio, if requested, is drawn as a diamond with horizontal limits at the
confidence limits and width inversely proportional to its standard
error.
}
\value{
  An object of class \code{meta.DSL} with \code{print}, \code{plot}, \code{funnelplot} and
  \code{summary} methods.
}
\author{Thomas Lumley}

\seealso{\code{\link{plot}},\code{\link{par}},\code{\link{meta.MH}},\code{\link{funnelplot}} }

\examples{
data(catheter)
b <- meta.DSL(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
              names=Name, subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
b
summary(b)
plot(b)
e <- meta.DSL(n.trt, n.ctrl, inf.trt, inf.ctrl, data=catheter,
              names=Name, subset=c(13,6,3,12,4,11,1,14,8,10,2))
e
summary(e)
##tasteless
plot(e, colors=meta.colors(summary="green",lines="purple",box="orange"))
}
\keyword{htest}
\keyword{hplot}
