\name{cv.dataSharing}
\alias{cv.dataSharing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross Validaate a Coupled Linear Model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv.dataSharing(predictors, response, groupsId, nfolds = 10, gamma = 2, tausGrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictors}{
%%     ~~Describe \code{predictors} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{groupsId}{
%%     ~~Describe \code{groupsId} here~~
}
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{tausGrid}{
%%     ~~Describe \code{tausGrid} here~~
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
function (predictors, response, groupsId, nfolds = 10, gamma = 2, 
    tausGrid) 
{
    stopifnot(nrow(predictors) == length(response))
    nTrain <- length(response)
    foldsId <- sample(rep(1:nfolds, length.out = nTrain))
    cvResults <- matrix(NA, nrow = nfolds, ncol = nrow(tausGrid))
    for (k in 1:nfolds) {
        testId <- which(foldsId == k)
        trainX <- predictors[-testId, ]
        trainY <- response[-testId]
        trainG <- groupsId[-testId]
        testX <- predictors[testId, ]
        testY <- response[testId]
        testG <- groupsId[testId]
        fittedModels <- apply(t(tausGrid), 2, function(myTaus) dataSharing(trainX, 
            trainY, trainG, gamma, myTaus))
        preditions <- apply(t(1:nrow(tausGrid)), 2, function(i) predict.dataSharing(fittedModels[, 
            i], testX, testG, tausGrid[i, ]))
        cvResults[k, ] <- apply(preditions, 2, function(yHat) mean((testY - 
            yHat)^2))
        print(paste("Fold", k, "best taus test MSE is", min(cvResults[k, 
            ])))
    }
    cv <- list()
    cv$cvm <- colMeans(cvResults)
    cv$cvsd <- apply(cvResults, 2, sd)
    cv$taus.min <- tausGrid[which.min(cv$cvm), ]
    cv$dataSharing.fit <- dataSharing(predictors, response, groupsId, 
        gamma, cv$taus.min)
    return(cv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
