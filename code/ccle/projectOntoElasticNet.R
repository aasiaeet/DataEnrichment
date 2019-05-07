######################################################################
# Data Sharing for Cancer Cell Line Encyclopedia (CCLE)
# Author: Amir Asiaee T. 
# Email: asiae002@umn.edu
######################################################################


evaluateElasticNet <- function(point, gamma){
  return(sum(abs(point)) + (gamma / 2) * sum(point^2)) 
}

projectOntoElasticNet <- function(projectMe, gamma, radius){
  if (evaluateElasticNet(projectMe, gamma) <= radius)
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
      if(s + deltaS - (rho + deltaRho) * evaluateElasticNet(currentValue, gamma) < radius * (1 + gamma * currentValue)**2){
        s <- s + deltaS
        rho <- deltaRho
        U <- belowIndexSet
      }
      else{
        U <- aboveIndexSet[-1]
      }
      sortedAbsProjectMe <- sortedAbsProjectMe[sortedAbsProjectMe$ix %in% U, ]
    }
    a <- gamma^2 * radius + gamma*rho / 2
    b <- 2 * gamma * radius + rho
    c <- radius - s 
    lambda <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
    return (sapply(projectMe, function(x) sign(x)* max(abs(x) - lambda,0) / (1 + lambda * gamma)))
  }
}

# for small gamma this implementation does not work? 
# print(projectOntoElasticNet(c(2,0,-2), .001, 1))


