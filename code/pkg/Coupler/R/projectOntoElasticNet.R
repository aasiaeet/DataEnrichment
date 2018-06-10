######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################
norm2sq <- function(x) (sum(x^2))

evaluateElasticNet <- function(point, gamma){
  return (norm(as.matrix(point), type = c("o")) + (gamma / 2) * norm2sq(as.matrix(point)))
}

projectOntoElasticNet <- function(projectMe, gamma, tau){
  if (evaluateElasticNet(projectMe, gamma) <= tau)
    return (projectMe)
  else{
    U <- 1:length(projectMe);
    s <- 0;
    rho <- 0
    
    absProjectMe <- abs(projectMe)
    sortedAbsProjectMe <- as.data.frame(sort(absProjectMe, index.return = TRUE))
    
    set.seed(123)
    while(length(U) != 0){
      currentIndex <- sample(1:length(U), 1)
      currentValue <- sortedAbsProjectMe$x[currentIndex]
      
      aboveIndexSet <- sortedAbsProjectMe$ix[currentIndex:length(sortedAbsProjectMe$ix)]
      aboveValueVec <- sortedAbsProjectMe$x[currentIndex:length(sortedAbsProjectMe$ix)]
      
      belowIndexSet <- sortedAbsProjectMe$ix[1:currentIndex - 1]
      belowValueVec <- sortedAbsProjectMe$x[1:currentIndex - 1]
      
      deltaRho <- length(aboveIndexSet)
      deltaS <- evaluateElasticNet(aboveValueVec, gamma)
      if(s + deltaS - (rho + deltaRho) * evaluateElasticNet(currentValue, gamma) < tau * (1 + gamma * currentValue)**2){
        s <- s + deltaS
        rho <- deltaRho
        U <- belowIndexSet
      }
      else{
        U <- aboveIndexSet[-1]
      }
      sortedAbsProjectMe <- sortedAbsProjectMe[sortedAbsProjectMe$ix %in% U, ]
    }
    a <- gamma**2 * tau + gamma*rho / 2
    b <- 2 * gamma * tau + rho
    c <- tau - s 
    lambda <- (-b + sqrt(b**2 - 4 * a * c)) / (2 * a)
    return (unlist(lapply(projectMe, function(x) sign(x)* max(abs(x) - lambda,0) / (1 + lambda * gamma))))
  }
}

print(projectOntoElasticNet(c(3,0,-3), .1, 1))
  

