#############################
# Projecting to l1 ball, R equivalent of the following MATLAB code by John Duchi:
# https://stanford.edu/~jduchi/projects/DuchiShSiCh08/ProjectOntoL1Ball.m
# Author: asiae002@umn.edu
#############################
projectOntoL1 <- function(projectMe, tau){
  stopifnot(tau >= 0)
  if(tau == 0){
    return(rep(0,length(projectMe)))
  }
  if(sum(abs(projectMe)) < tau)
    return(projectMe)
  u <- sort(abs(projectMe), decreasing = TRUE)
  sv <- cumsum(u)
  rho <- max(which( u > (sv - tau) / 1:length(u) ))
  theta <- max(0, (sv[rho] - tau) / rho)
  return( sign(projectMe) * pmax( abs(projectMe) - theta, rep(0, length(projectMe)) ))
}

# print(projectOntoL1(1:3,0))
