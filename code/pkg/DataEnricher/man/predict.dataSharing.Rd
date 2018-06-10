\name{predict.dataSharing}
\alias{predict.dataSharing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Predictions from a Coupled Linear Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.dataSharing(fittedModel, x, groups, taus, addIntercept = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fittedModel}{
%%     ~~Describe \code{fittedModel} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{taus}{
%%     ~~Describe \code{taus} here~~
}
  \item{addIntercept}{
%%     ~~Describe \code{addIntercept} here~~
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
function (fittedModel, x, groups, taus, addIntercept = TRUE) 
{
    if (addIntercept == TRUE) {
        x <- cbind(x, rep(1, nrow(x)))
    }
    learnedBeta <- matrix(fittedModel, ncol = 3)
    predictedY <- rep(0, nrow(x))
    for (g in 1:2) {
        groupIndex <- groups == g
        predictedY[groupIndex] <- x[groupIndex, ] \%*\% (learnedBeta[, 
            g] + learnedBeta[, 3])
    }
    return(predictedY)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
