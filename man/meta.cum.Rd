\name{meta.cum}
\alias{meta.cum}
\alias{print.meta.cum}
\alias{plot.meta.cum}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cumulative meta-analysis of binary data}
\description{
A cumulative meta-analysis plot shows how evidence has accumulated over
time.  The \code{i}th line in the cumulative meta-analysis plot is the
summary produced by a meta-analysis of the first \code{i} trials.
}
\usage{
meta.cum(ntrt, nctrl, ptrt, pctrl, conf.level = 0.95, names = NULL, data= NULL, subset = NULL,
 na.action =get(getOption("na.action")), method = c("meta.MH", "meta.DSL"), statistic = "OR")
plot.meta.cum(x, conf.level = NULL, colors = meta.colors(), xlab = NULL, 
    summary.line = TRUE, summary.conf = FALSE, ...) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{ntrt}{Number of subjects in treated/exposed group}
 \item{nctrl}{Number of subjects in control group }
 \item{ptrt}{Number of events in treated/exposed group}
 \item{pctrl}{Number of events in control group}
 \item{names}{names or labels for studies}
 \item{data}{data frame to interpret variables}
 \item{subset}{subset of studies to include}
 \item{na.action}{How to handle missing values}
  \item{method}{Which meta-analysis method to use}
  \item{statistic}{"OR" for odds ratio or "RR" for relative risk. }
 \item{x}{a \code{meta.cum} object}
 \item{\dots}{other graphical arguments}
 \item{conf.level}{Coverage for confidence intervals}
 \item{colors}{see \code{\link{meta.colors}}}
 \item{xlab}{X-axis label}
 \item{summary.line}{Plot a vertical line at the final summary value?}
 \item{summary.conf}{Plot vertical lines at the final confidence
   interval limits?}
}


\value{
  Object of class \code{meta.cum}.
}
\seealso{\code{\link{meta.MH}},\code{\link{meta.DSL}}}
\examples{
data(cochrane)
steroid<-meta.cum(n.trt,n.ctrl,ev.trt,ev.ctrl,names=name,data=cochrane,statistic="RR",method="meta.MH")
plot(steroid)
}
\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
