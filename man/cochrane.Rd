\name{cochrane}
\alias{cochrane}
\non_function{}
\title{Data for Cochrane Collaboration logo }
\usage{data(cochrane)}
\description{
Data from randomised trials before 1980 of corticosteroid therapy in premature
labour and its effect on neonatal death. 
}
\format{
This data frame contains the following columns:
\describe{
\item{name}{Identifier for the study
}
\item{ev.trt}{
Number of deaths in the treated group
}
\item{n.trt}{
Number in the treated group
}
\item{ev.ctrl}{
Number of deaths in the control group
}
\item{n.ctrl}{
Number in the control group
}
}
}
\details{
This meta-analysis, if done, would likely have resulted in the treatment being
widely used a decade earlier than it was, saving many lives. The graph
is part of the logo of the Cochrane Collaboration, a group aiming to
perform systematic reviews of the entire clinical trial literature.
} \source{
\url{http://www.cochrane.org}
}
\examples{
data(cochrane)
steroid <- meta.MH(n.trt, n.ctrl, ev.trt, ev.ctrl,
                   names=name, data=cochrane)
plot(steroid, col=meta.colors("RoyalBlue"))
}
\keyword{datasets}
