\name{catheter}
\title{Meta-analysis of antibacterial catheter coating }
\usage{data(catheter)}
\alias{catheter}
\description{Data on the effectiveness of silver sulfadiazine coating on
  venous catheters for preventing bacterial colonisation of the catheter
  and bloodstream infection}
\source{Veenstra D et al (1998) "Efficacy of Antiseptic Impregnated
Central Venous Catheters in Preventing Nosocomial Infections: A
Meta-analysis" JAMA 281:261-267}
\format{A data.frame with 8 variables giving information about 16
  controlled trials of antibacterial-coated venous catheters
  \tabular{rl}{
    Name      :  \tab Name of principal author\cr
    n.trt     :  \tab number of coated catheters\cr
    n.ctrl    : \tab number of standard catheters\cr
    col.trt   : \tab number of coated catheters colonised\cr
    \tab by bacteria\cr
    col.ctrl  :\tab number of standard catheters colonised\cr
    \tab by bacteria\cr
    inf.trt   : \tab number of coated catheters resulting in \cr
    \tab bloodstream infection\cr
    inf.ctrl  :\tab number of standard catheters resulting in\cr
    \tab bloodstream infection\cr
    or       :\tab Odds ratio\cr
  }
}
\examples{
library(rmeta)
data(catheter)
a <- meta.MH(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
             names=Name, subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
b <- meta.DSL(n.trt, n.ctrl, col.trt, col.ctrl, data=catheter,
              names=Name,subset=c(13,6,5,3,7,12,4,11,1,8,10,2))
a
b
summary(a)
summary(b)
plot(a)
plot(b)
}
\keyword{datasets}
