\name{meta.MH}
\title{Fixed effects (Mantel-Haenszel) meta-analysis}
\usage{
meta.MH(ntrt, nctrl, ptrt, pctrl, names=NULL, data=NULL, subset=NULL)
summary.meta.MH(object)
plot.meta.MH(object,summary=T,summlabel="Summary",\dots)
}
\alias{meta.MH}
\alias{summary.meta.MH}
\alias{plot.meta.MH}
%- Also NEED an `\alias' for EACH other function documented here.
\arguments{
 \item{ntrt}{Number of subjects in treated/exposed group}
 \item{nctrl}{Number of subjects in control group }
 \item{ptrt}{Number of events in treated/exposed group}
 \item{pctrl}{Number of events in control group}
 \item{names}{names or labels for studies}
 \item{data}{data frame to interpret variables}
 \item{subset}{subset of studies to include}
 \item{object}{a \code{meta.MH} object}
 \item{summary}{Plot the summary odds ratio?}
 \item{summlabel}{Label for the summary odds ratio}
 \item{\dots}{other graphical arguments}
}
\description{
Computes the individual odds ratios, the Mantel-Haenszel summary odds
ratio, and Woolf's test for heterogeneity. The \code{print} method gives
the summary and test for heterogeneity; the \code{summary} method also
gives all the individual odds ratios and confidence intervals.

The \code{plot} method draws a standard meta-analysis plot. The
confidence interval for each study is given by a horizontal line, and
the point estimate is given by a square whose height is inversely
proportional to the standard error of the estimate.  The summary odds
ratio, if requested, is drawn as a diamond with horizontal limits at the
confidence limits and width inversely proportional to its standard
error.
}
\value{
  An object of class \code{meta.MH} with \code{print}, \code{plot} and
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
}
\references{}
\author{Thomas Lumley}
\note{ There are at least two other ways to do a fixed effects
  meta-analysis of binary data. Peto's method is a computationally
  simpler approximation to the Mantel-Haenszel approach. It is also
  possible to weight the individual odds ratios according to their
  estimated variances.  The Mantel-Haenszel method is superior if there
  are trials with small numbers of events (less than 5 or so in either group)}

\seealso{\code{\link{plot}},\code{\link{par}},\code{\link{meta.DSL}} }

\examples{
data(catheter)
a<-meta.MH(n.trt,n.ctrl,col.trt,col.ctrl,data=catheter,names=Name,subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
a
summary(a)
plot(a)
d<-meta.MH(n.trt,n.ctrl,inf.trt,inf.ctrl,data=catheter,names=Name,subset=c(13,6,3,12,4,11,1,14,8,10,2))
d
summary(d)
plot(d)

}
\keyword{htest}
\keyword{hplot}
