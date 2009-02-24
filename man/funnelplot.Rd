\name{funnelplot}
\alias{funnelplot}
\alias{funnelplot.meta.MH}
\alias{funnelplot.meta.DSL}
\alias{funnelplot.meta.summaries}
\alias{funnelplot.default}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Funnel plot for publication bias}
\description{
  Plots the treatment difference for trials against the size of the
  trial (or other specified variable).  Asymmetry in the plot often
  indicates publication bias.  Generic, with methods for meta-analysis
  objects.
}
\usage{
funnelplot(x,\dots)
\method{funnelplot}{default}(x, se, size=1/se, summ=NULL,
        xlab="Effect", ylab="Size", colors=meta.colors(),
        conf.level=0.95, plot.conf=FALSE,
        zero=NULL, mirror=FALSE, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Treatment difference}
  \item{se}{Standard error of \code{x}}
  \item{size}{Variable for the vertical axis}
  \item{summ}{summary treatment difference}
  \item{xlab}{x-axis label }
  \item{ylab}{y-axis label}
  \item{colors}{list of colors for components of the plot}
  \item{conf.level}{For confidence interval plotting }
  \item{plot.conf}{Plot confidence intervals instead of just points?}
  \item{zero}{location of a null hypothesis line}
  \item{mirror}{Add points reflected around \code{summ}?}
  \item{\dots}{further arguments to be passed from or to methods.}
}
\details{
  With the default value of \code{size} the plot should appear as a
  upwards-pointing funnel shape.  Publication bias often causes one side
  of the funnel to be trimmed near the base.  The \code{mirror} plot
  creates a symmetric funnel by reflecting the plot around the
  \code{summ} value.  In the presence of publication bias the added
  points will separate from the real studies.  
}
\value{
  Used for its side-effect.
}
\author{Thomas Lumley}

\seealso{
  \code{\link{meta.DSL}},
  \code{\link{meta.MH}},
  \code{\link{meta.summaries}},
  \code{\link{metaplot}}
}
\examples{
data(catheter)
a <- meta.MH(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
             names=Name, subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
funnelplot(a$logOR, a$selogOR)
funnelplot(a$logOR, a$selogOR,
           plot.conf=TRUE, summ=a$logMH, mirror=TRUE)
funnelplot(a, plot.conf=TRUE)
}
\keyword{hplot}
