

stratify <- function(dat, groups, epsilon=0.000001, maxiter=1000, try=20){
  # dat - must samples in columns and edges in rows ..
  # groups - the number of groups to stratify to
  # epsilon - when to stop the EM algorithm
  # maxiter - when to stop the EM algorithm
  # try - number of tries if the initial state fails to converge at maxiters
  stopifnot(class(dat) == "matrix")
  stopifnot(class(groups) == "numeric", length(groups) == 1)
  stopifnot(class(epsilon) == "numeric")
  stopifnot(class(maxiter) == "numeric")
  stopifnot(class(try) == "numeric")
  res0 <- 0
  while(class(res0) != "list" & try > 0) {
    res0 <- .Call("stratem", dat, groups, epsilon, maxiter, PACKAGE = "hyperstrat" )
    try <- try-1
  }
  if(class(res0) != "list") {
    print("\n****************\nWARNING FAILED STRATIFICATION\n****************\n")
  }
  res0
}


printworkbar <- function() {
  # prints a status bar showing the progress of the stratification
  # if doing multiple strats.
  cat("\nStratifying...\n||")
  for (i in 1:50) cat ("-")
  cat("||\n  ")
}



multistrat <- function(dat, groups, times, collapse=T,
                       epsilon=0.000001, maxiter=1000, try=20) {
  # dat - must have samples in columns and edges in rows ..
  # groups - the number of groups to stratify to
  # times - the number of stratification replicates to do
  # collapse - merge the strat. replicates? or return a list.
  # epsilon - when to stop the EM algorithm
  # maxiter - when to stop the EM algorithm
  # try - number of tries if the initial state fails to converge at maxiters
  stopifnot(class(dat) == "matrix")
  stopifnot(class(groups) == "numeric", length(groups) == 1)
  stopifnot(class(times) == "numeric", length(times) == 1)
  printworkbar()
  reslist <- vector("list", times)
  if (times < 50) {chunksize <- as.integer(50/times); printfreq <- 1}
  else {chunksize <- 1; printfreq <- as.integer(50/times)+1;} # when to print a progress update
  for (i in 1:times) {
    # produce a stratification
    reslist[[i]] <- stratify(dat,groups,epsilon,maxiter,try)
    if (i %% printfreq == 0) {for (k in 1:chunksize) {cat("=")}}
  }
  cat(">\n  Stratifications Done!\n")
  if (collapse) {
    cat(">\n  Collapsing ... \n")
    reslist <- reslist[sapply(reslist, function(b) !is.null(b))]
    res1 <- collapsestrats(reslist)
    return(res1)
  } else {
    return(reslist)
  }
}



parmultistrat <- function(dat, groups, times=1, threads=1, collapse=F,
                          epsilon=0.000001, maxiter=1000, try=20) {
  # dat - must have samples in columns and edges in rows ..
  # groups - the number of groups to stratify to
  # times - the number of stratification replicates to do
  # threads - number of threads to perform in parallel .. uses multicore
  # collapse - merge the strat. replicates? or return a list.
  # epsilon - when to stop the EM algorithm
  # maxiter - when to stop the EM algorithm
  # try - number of tries if the initial state fails to converge at maxiters
  stopifnot(class(dat) == "matrix")
  stopifnot(class(groups) == "numeric", length(groups) == 1)
  stopifnot(class(times) == "numeric", length(times) == 1)
  require(multicore)
  print("Starting parallel processes...")
  reslist <- mclapply(X=1:times, function(i) {
    cat(i,"-:-"); stratify(dat,groups,epsilon,maxiter,try);
  }, mc.cores=threads)
  print("Done!")
  if (collapse) {
    cat(">\n  Collapsing ... \n")
    reslist <- reslist[sapply(reslist, function(b) !is.null(b))]
    res1 <- collapsestrats(reslist)
    return(res1)
  } else {
    return(reslist)
  }
}

