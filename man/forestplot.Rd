\name{forestplot}
\alias{forestplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Forest plots }
\description{
Draw a forest plot together with a table of text.
}
\usage{
forestplot(labeltext, mean, lower, upper, align = NULL, is.summary =
FALSE, clip = c(-Inf, Inf), xlab = "", zero = 0, graphwidth = unit(2,
"inches"), col = meta.colors(), xlog = FALSE, xticks=NULL, boxsize=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labeltext}{Matrix of strings or \code{NA}s for blank spaces }
  \item{mean}{Vector of centers of confidence intervals (or \code{NA}s for
    blank space) }
  \item{lower}{Vector of lower ends of confidence intervals }
  \item{upper}{Vector of upper ends of confidence intervals  }
  \item{align}{Vector giving alignment (\code{l},\code{r},\code{c}) for
    columns of table }
  \item{is.summary}{Vector of logicals. Summary lines have bold text and
    diamond confidence intervals.
  }
  \item{clip}{Lower and upper limits for clipping confidence intervals
    to arrows }
  \item{xlab}{x-axis label }
  \item{zero}{x-axis coordinate for zero line }
  \item{graphwidth}{Width of confidence interval graph }
  \item{col}{See \code{\link{meta.colors}}}
  \item{xlog}{If \code{TRUE}, x-axis tick marks are exponentiated}
  \item{xticks}{Optional user-specified x-axis tick marks. Specify \code{NULL}
    to use the defaults, \code{numeric(0)} to omit the x-axis.}
  \item{boxsize}{Override the default box size based on precision}
  \item{\dots}{Not used.}
}
\details{
This function is more flexible than \code{\link{metaplot}} and the
\code{plot} methods for meta-analysis objects, but requires more work by
the user.

In particular, it allows for a table of text, and clips confidence
intervals to arrows when they exceed specified limits.
}
\value{
None
}

\seealso{\code{metaplot}}
\examples{
data(cochrane)
steroid <- meta.MH(n.trt, n.ctrl, ev.trt, ev.ctrl,
                        names=name, data=cochrane)

tabletext<-cbind(c("","Study",steroid$names,NA,"Summary"),
                 c("Deaths","(steroid)",cochrane$ev.trt,NA,NA),
                 c("Deaths","(placebo)", cochrane$ev.ctrl, NA,NA),
                 c("","OR",format(exp(steroid$logOR),digits=2),NA,format(exp(steroid$logMH),digits=2))
                 )

m<- c(NA,NA,steroid$logOR,NA,steroid$logMH)
l<- m-c(NA,NA,steroid$selogOR,NA,steroid$selogMH)*2
u<- m+c(NA,NA,steroid$selogOR,NA,steroid$selogMH)*2
forestplot(tabletext,m,l,u,zero=0,is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
   clip=c(log(0.1),log(2.5)), xlog=TRUE,
   col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))

forestplot(tabletext,m,l,u,zero=0,is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
   clip=c(log(0.1),log(2.5)), xlog=TRUE, boxsize=0.75,
   col=meta.colors(box="royalblue",line="darkblue", summary="royalblue"))

}
\keyword{hplot}% at least one, from doc/KEYWORDS

