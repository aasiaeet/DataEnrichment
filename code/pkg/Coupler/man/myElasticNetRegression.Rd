\name{myElasticNetRegression}
\alias{myElasticNetRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
I Still Don't KNowe
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
myElasticNetRegression(data, givenLambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{givenLambda}{
%%     ~~Describe \code{givenLambda} here~~
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
function (data, givenLambda) 
{
    response <- as.matrix(data[, "y", drop = FALSE])
    predictors <- as.matrix(data[, !(names(data) \%in\% c("y", 
        "colors"))])
    best.predictors <- predictors
    print("Starting Elastic Net Regression ...")
    if (missing(givenLambda)) {
        print("There is no Lambda! Let's cross-validate!")
        tmpCv1005 <- cv.glmnet(best.predictors, response, nfolds = 10, 
            alpha = 0.5)
        print("Done with the cross-validation of elastic net!")
    }
    else {
        print("There is a Lambda! Let's use it!")
        tmpCv1005 <- glmnet(best.predictors, response, alpha = 0.5, 
            lambda = givenLambda)
        print("Done with the elastic net!")
    }
    return(tmpCv1005)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
