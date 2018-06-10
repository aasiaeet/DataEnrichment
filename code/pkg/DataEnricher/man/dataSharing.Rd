\name{dataSharing}
\alias{dataSharing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Coupled Linear Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dataSharing(x, y, groups, gamma, taus)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{taus}{
%%     ~~Describe \code{taus} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, groups, gamma, taus) 
{
    eta <- 1/(10 * nrow(x))
    x <- cbind(x, rep(1, nrow(x)))
    BetaPrv <- matrix(0, ncol(x), 3)
    marginalImp <- rep(0, length(y))
    BetaNex <- BetaPrv
    meanSqErrOld <- 0
    for (i in 1:nStepsGD) {
        eta <- (1/(1 + i)) * eta
        gradient <- 0
        for (g in 1:2) {
            groupIndex <- groups == g
            scaledGradient <- sqrt(eta * (1/sum(groupIndex))) * 
                t(x[groupIndex, ]) \%*\% (y[groupIndex] - x[groupIndex, 
                  ] \%*\% (BetaPrv[, g] + BetaPrv[, 3]))
            BetaNex[, g] <- BetaPrv[, g] + scaledGradient
            gradient <- gradient + scaledGradient
            BetaNex[, g] <- projectOntoElasticNet(BetaNex[, g], 
                gamma, as.double(taus[g]))
        }
        BetaNex[, 3] <- BetaPrv[, 3] + gradient
        BetaNex[, 3] <- projectOntoElasticNet(BetaNex[, 3], gamma, 
            as.double(taus[3]))
        predictedY <- predict.dataSharing(BetaNex, x, groups, 
            taus, addIntercept = FALSE)
        meanSqErrNew <- mean((y - predictedY)^2)
        if (abs(meanSqErrNew - meanSqErrOld) < stoppingCriteria) {
            break
        }
        meanSqErrOld <- meanSqErrNew
        BetaPrv <- BetaNex
    }
    return(BetaNex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
