\name{meta.DSL}
\title{Random effects (DerSimonian-Laird) meta-analysis}
\usage{
meta.DSL(ntrt, nctrl, ptrt, pctrl, names=NULL, data=NULL,
 subset=NULL,conf.level=0.95)
summary.meta.DSL(object,conf.level=NULL)
plot.meta.DSL(object,summary=T,summlabel="Summary",\dots,conf.level=NULL,
colors=meta.colors())
}
\alias{meta.DSL}
\alias{summary.meta.DSL}
\alias{plot.meta.DSL}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ntrt}{Number of subjects in treated/exposed group}
 \item{nctrl}{Number of subjects in control group }
 \item{ptrt}{Number of events in treated/exposed group}
 \item{pctrl}{Number of events in control group}
 \item{names}{names or labels for studies}
 \item{data}{data frame to interpret variables}
 \item{subset}{subset of studies to include}
 \item{object}{a \code{meta.DSL} object}
 \item{summary}{Plot the summary odds ratio?}
 \item{summlabel}{Label for the summary odds ratio}
 \item{\dots}{other graphical arguments}
 \item{conf.level}{Coverage for confidence intervals}
 \item{colors}{see \code{\link{meta.colors}}}
}
\description{
Computes the individual odds ratios, the summary odds ratio, the random
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
  An object of class \code{meta.DSL} with \code{print}, \code{plot} and
  \code{summary} methods and components:
  \item{logOR}{log odds ratios for individual studies}
  \item{selogOR}{standard errors for log odds ratios}
  \item{logMH}{log of Mantel-Haenszel summary odds ratio}
  \item{selogMH}{standard of summary log odds ratio}
  \item{MHtest}{Mantel-Haenszel chisquare and p-value testing the
    hypothesis that the summary odds ratio is 1}
  \item{het}{Woolf's chisquare for heterogeneity, its degrees of freedom
    and p-value}
  \item{call}{A copy of the function call}
  \item{names}{A copy of the vector of names}
  \item{omitted}{logical vector: was the study omitted because of 0/Inf
    odds ratio}
  \item{tau2}{Estimated random effects variance}
}
\author{Thomas Lumley}

\seealso{\code{\link{plot}},\code{\link{par}},\code{\link{meta.MH}} }

\examples{
data(catheter)
b<-meta.DSL(n.trt,n.ctrl,col.trt,col.ctrl,data=catheter,names=Name,subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
b
summary(b)
plot(b)
e<-meta.DSL(n.trt,n.ctrl,inf.trt,inf.ctrl,data=catheter,names=Name,subset=c(13,6,3,12,4,11,1,14,8,10,2))
e
summary(e)
#tasteless
plot(e,colors=meta.colors(summary="green",lines="purple",box="orange"))
}
\keyword{htest}
\keyword{hplot}
