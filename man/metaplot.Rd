\name{metaplot}
\alias{metaplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Meta-analysis plot (forest plot)}
\description{
  Plot confidence intervals with boxes indicating the sample
  size/precision and optionally a diamond indicating a summary
  confidence interval. This function is usually called by \code{plot}
  methods for meta-analysis objects.
}
\usage{
metaplot(mn, se, nn=NULL, labels=NULL, conf.level=0.95,
         xlab="Odds ratio", ylab="Study Reference",xlim=NULL,
         summn=NULL, sumse=NULL, sumnn=NULL, summlabel="Summary",
         logeffect=FALSE, lwd=2, boxsize=1,
         zero=as.numeric(logeffect), colors=meta.colors(),
         xaxt="s", logticks=TRUE,  ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mn}{point estimates from studies}
  \item{se}{standard errors of \code{mn}}
  \item{nn}{precision: box ares is proportional to this. \code{1/se^2}
    is the default}
  \item{labels}{labels for each interval}
  \item{conf.level}{Confidence level for confidence intervals}
  \item{xlab}{label for the point estimate axis}
  \item{ylab}{label for the axis indexing the different studies}
  \item{xlim}{the range for the x axis.}
  \item{summn}{summary estimate}
  \item{sumse}{standard error of summary estimate}
  \item{sumnn}{precision of summary estimate}
  \item{summlabel}{label for summary estimate }
  \item{logeffect}{\code{TRUE} to display on a log scale }
  \item{lwd}{line width}
  \item{boxsize}{Scale factor for box size}
  \item{zero}{"Null" effect value}
  \item{xaxt}{use \code{"n"} for no x-axis (to add a customised one)}
  \item{logticks}{if \code{TRUE} and \code{logscale}, have tick values
    approximately equally spaced on a log scale}.
  \item{colors}{see \code{\link{meta.colors}}}
  \item{\dots}{Other graphical parameters }
}

\value{
  This function is used for its side-effect.
}

\seealso{
  \code{\link{forestplot}} for more flexible plots
  
  \code{\link{plot.meta.DSL}},
  \code{\link{plot.meta.MH}},
  \code{\link{plot.meta.summaries}}
}

\examples{
data(catheter)
a <- meta.MH(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
             names=Name, subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
metaplot(a$logOR, a$selogOR, nn=a$selogOR^-2, a$names,
         summn=a$logMH, sumse=a$selogMH, sumnn=a$selogMH^-2,
         logeffect=TRUE)
metaplot(a$logOR, a$selogOR, nn=a$selogOR^-2, a$names,
         summn=a$logMH, sumse=a$selogMH, sumnn=a$selogMH^-2,
         logeffect=TRUE,logticks=FALSE)
## angry fruit salad
metaplot(a$logOR, a$selogOR, nn=a$selogOR^-2, a$names,
         summn=a$logMH, sumse=a$selogMH, sumnn=a$selogMH^-2,
         logeffect=TRUE, colors=meta.colors(box="magenta",
             lines="blue", zero="red", summary="orange",
             text="forestgreen"))
}
\keyword{hplot}

