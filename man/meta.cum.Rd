\name{cummeta}
\alias{cummeta}
\alias{cummeta.summaries}
\alias{print.meta.cum}
\alias{plot.meta.cum}
\alias{print.summary.meta.cum}
\alias{summary.meta.cum}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Cumulative meta-analysis of binary data}
\description{
A cumulative meta-analysis plot shows how evidence has accumulated over
time.  The \code{i}th line in the cumulative meta-analysis plot is the
summary produced by a meta-analysis of the first \code{i} trials.
}
\usage{
cummeta(ntrt, nctrl, ptrt, pctrl, conf.level = 0.95, names = NULL, data=
NULL,  subset = NULL, na.action = na.fail, method = c("meta.MH", "meta.DSL"), statistic = "OR")
cummeta.summaries(effects,stderrs, conf.level = 0.95, names = NULL,weights=NULL, data= NULL, subset = NULL, na.action = get(getOption("na.action")), method = c("fixed", "random"), logscale=TRUE)
\method{plot}{meta.cum}(x, conf.level = NULL, colors = meta.colors(), xlab = NULL, 
    summary.line = TRUE, summary.conf = FALSE, main="Cumulative meta-analysis", lwd=1, ...)
\method{summary}{meta.cum}(object ,conf.level=NULL,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{ntrt}{Number of subjects in treated/exposed group}
 \item{nctrl}{Number of subjects in control group }
 \item{ptrt}{Number of events in treated/exposed group}
 \item{pctrl}{Number of events in control group}
 \item{effects}{Difference between control and treatment group}
 \item{stderrs}{Standard errors of \code{effects}}
 \item{weights}{Study weights (see \code{\link{meta.summaries}})}
 \item{names}{names or labels for studies} 
 \item{data}{data frame to interpret variables}
 \item{subset}{subset of studies to include}
 \item{na.action}{How to handle missing values}
  \item{method}{Which meta-analysis method to use}
  \item{statistic}{"OR" for odds ratio or "RR" for relative risk. }
  \item{logscale}{The \code{effects} values are the logarithms of actual
    effects (for plotting)}
  \item{x, object}{a \code{meta.cum} object}
  \item{\dots}{other graphical arguments for \code{\link{metaplot}}}
  \item{conf.level}{Coverage for confidence intervals}
  \item{colors}{see \code{\link{meta.colors}}}
  \item{xlab}{X-axis label}
  \item{summary.line}{Plot a vertical line at the final summary value?}
  \item{summary.conf}{Plot vertical lines at the final confidence
    interval limits?}
  \item{main,lwd}{graphical parameters}
}


\value{
  Object of class \code{meta.cum}.
}
\seealso{\code{\link{meta.MH}},\code{\link{meta.DSL}}}
\examples{
data(cochrane)
steroid<-cummeta(n.trt,n.ctrl,ev.trt,ev.ctrl,names=name,data=cochrane,statistic="RR",method="meta.MH")
plot(steroid)
summary(steroid)
data(catheter)
b <- meta.DSL(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
              names=Name, subset=c(13,6,5,3,12,4,11,1,8,10,2))
d <- cummeta.summaries(b$logs, b$selogs, names=b$names,
                    method="random", logscale=TRUE)
plot(d,summary.conf=TRUE)
summary(d)
}
\keyword{htest}% at least one, from doc/KEYWORDS
\keyword{hplot}% __ONLY ONE__ keyword per line
