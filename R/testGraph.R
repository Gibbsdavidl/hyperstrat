

generate <- function(params=16, nodes=60, groups=4, paramsPerGroup=5) {

  # generate the parameters
  theta <- matrix(nrow=groups, ncol=params, data=0)
  for (i in 1:groups) {
    theta[i, sample(1:params, size=paramsPerGroup, replace=F)] <- runif(paramsPerGroup)
  }
  # make the group assignments
  q <- matrix(nrow=nodes, ncol=groups, data=0)
  nodegrps <- sample(1:groups, size=nodes, replace=T)
  for (i in 1:groups) {
    for (j in 1:nodes) {
      if (nodegrps[j] == i) {
        q[j,i] <- 1
      }
    }
  }
  # make the adjacency matrix
  a <- matrix(ncol=params, nrow=nodes, data=0)
  for (i in 1:nodes) {
    for (j in 1:params) {
      if (runif(1) < theta[nodegrps[i], j]) {
        a[i,j] <- 1
      }
    }
  }
  list(a, theta, q, nodegrps)
}


teststrat <- function(g=NULL,params=16, nodes=60, groups=4, paramsPerGroup=5,
                      single=T, threads=1, times=10) {
  f <- function(a) which(a == max(a))

  if (is.null(g)) {
    g <- generate(params, nodes, groups, paramsPerGroup)
  }
  if (single) {
    s <- stratify(t(g[[1]]), groups)
  } else {
    s <- parmultistrat(t(g[[1]]),groups, times, threads, collapse=T)
  }
  i <- optimizeColOrder(g[[3]], s[[1]])
  g[[1]] <- g[[1]][,i[[1]]]; g[[2]] <- g[[2]][i[[1]],]; g[[3]] <- g[[3]][,i[[1]]];
  s[[1]] <- s[[1]][,i[[2]]]; s[[2]] <- s[[2]][i[[2]],];
  res0 <- comparestrats(s, g)
  err <- sum((s[[1]] - g[[3]])^2)
  g[[4]] <- apply(g[[3]],1,f)
  estgrps <- apply(s[[1]],1,f)
  esterr <- sum(estgrps != g[[4]])

  print("Actual Groups")
  print(g[[4]])
  print("Estimated Groups")
  print(estgrps)
  print("Number Group Estimates in Error")
  print(esterr)
  print("Error in Group Assignments")
  print(err)
  list(err,res0,s,g,estgrps)
}




  


#stratify <- function(dat, groups){
#optimizeColOrder <- function(table1, table2) { 
#multistrat <- function(dat, groups, times, collapse=T) {
#collapsematrices <- function(reslist) {  
#collapsestrats <- function(reslist) {  
#parmultistrat <- function(dat, groups, times=1, threads=1, collapse=F) {
#comparestrats <- function(strat1, strat2) {
#findGroupN <- function(dat, groupRange, times, threads=1) {
#pkl <- function(a,b,k,l,m) {
#pk <- function(a,k,m) {
#mi <- function(strat1, strat2) {
