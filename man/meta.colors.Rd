\name{meta.colors}
\alias{meta.colors}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Control colours in meta-analysis plot}
\description{
 Wrapper function for specifying colours to meta-analysis plots
}
\usage{
meta.colors(all.elements, box="black", lines="gray", summary="black", zero="lightgray", mirror="lightblue",text="black")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{all.elements}{if present, overrides other arguments}
  \item{box}{Colour of sample size box }
  \item{lines}{Colour of confidence intervals }
  \item{summary}{Colour of summary estimate}
  \item{zero}{Colour of null hypothesis line}
  \item{mirror}{Colour of reflected points (in \code{\link{funnelplot}})
  }
  \item{text}{Colour of labels}
}

\value{
 a list of colors
}

\seealso{\code{\link{plot.meta.MH}},\code{\link{plot.meta.DSL}},\code{\link{plot.meta.summaries}},\code{\link{funnelplot}},\code{\link{metaplot}}}

\examples{

}
\keyword{hplot }%-- one or more ...
