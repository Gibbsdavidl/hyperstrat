
comparestrats <- function(strat1, strat2) {
  stopifnot(class(strat1) == "list") 
  stopifnot(class(strat2) == "list")
  q1 <- strat1[[1]]
  q2 <- strat2[[1]]
  stopifnot(class(q1) == "matrix")
  stopifnot(class(q2) == "matrix")
  x <- .Call("compareStrats", q1, q2, PACKAGE = "hyperstrat")
  x
}


findGroupN <- function(dat, groupRange, times, threads=1) {
  # dat - has samples in columns and edges in rows ..
  # groupRange - vector of groups to try
  # times - number of times to stratify per group
  # threads - number of parallel threads
  stopifnot(class(groupRange) == "numeric" | class(groupRange) == "integer")  
  res0 <- vector("list", length(groupRange))
  i <- 1
  for (g in groupRange) {
    print("**********************")
    cat  ("     Group ", g, "  \n")
    print("**********************")
    res0[[i]] <- parmultistrat(dat, g, times, threads, collapse=T, epsilon=0.00000001) 
    i <- i+1
  }
  mat <- mat.or.vec(length(groupRange),length(groupRange))
  for (i in 1:length(groupRange)) {
    for (j in 1:length(groupRange)) {
      mat[i,j] <- mi(res0[[i]], res0[[j]])[[1]]
    }
  }
  #diag(mat) <- 0
  list(mat, apply(mat, 2, mean))
}


pkl <- function(a,b,k,l,m) {
  sum(a[,k]*b[,l])/m
}

pk <- function(a,k,m) {
  sum(a[,k])/m
}

mi <- function(strat1, strat2) {
  stopifnot(class(strat1) == "list") 
  stopifnot(class(strat2) == "list")
  q1 <- strat1[[1]]
  q2 <- strat2[[1]]
  stopifnot(class(q1) == "matrix")
  stopifnot(class(q2) == "matrix")
  ni <- ncol(q1)
  nj <- ncol(q2)
  top <- 0
  bot1 <- 0
  bot2 <- 0
  jent <- 0
  m <- nrow(q1)

  jp <- 0

  for (k in 1:ni) {
    for (l in 1:nj) {
      if(pkl(q1,q2,k,l,m) > 0) { 
        jp <- jp + pkl(q1,q2,k,l,m)
        top = top + pkl(q1,q2,k,l,m)*log2(pkl(q1,q2,k,l,m)/(pk(q1,k,m)*pk(q2,l,m)))
      }
    }
  }

  for (k in 1:ni) {
    for (l in 1:nj) {
      if(pkl(q1,q2,k,l,m) > 0) { 
        jent = jent + -1 * pkl(q1,q2,k,l,m)*log2(pkl(q1,q2,k,l,m))
      }
    }
  }

  
  for (k in 1:ni) {
    if (pk(q1,k,m) > 0) {
      bot1 <- bot1 + -1*pk(q1,k,m)*log2(pk(q1,k,m))
    }
  }

  for (l in 1:nj) {
    if (pk(q2,l,m) > 0) {
      bot2 <- bot2 + -1*pk(q2,l,m)*log2(pk(q2,l,m))
    }
  }
  
  list(NormMI=top/(min(bot1,bot2)), MI=(bot1+bot2-jent), MI2=top, Jentro=jent, H1=bot1, H2=bot2, JP=jp)
}

exprmi <- function(strat1, strat2) {
  stopifnot(class(strat1) == "list") 
  stopifnot(class(strat2) == "list")
  q1 <- strat1[[1]]
  q2 <- strat2[[1]]
  stopifnot(class(q1) == "matrix")
  stopifnot(class(q2) == "matrix")
  ni <- ncol(q1)
  nj <- ncol(q2)
  Sij <- 0
  Sii <- 
  Sjj <- 0
  m <- nrow(q1)

  for (k in 1:ni) {
    for (l in 1:nj) {
      if(pkl(q1,q2,k,l,m) > 0) { 
        Sij = Sij + pkl(q1,q2,k,l,m)*log2(pkl(q1,q2,k,l,m)/(pk(q1,k,m)*pk(q2,l,m)))
      }
    }
  }

  for (k in 1:ni) {
    for (l in 1:ni) {
      if(pkl(q1,q1,k,l,m) > 0) { 
        Sii = Sii + pkl(q1,q1,k,l,m)*log2(pkl(q1,q1,k,l,m)/(pk(q1,k,m)*pk(q1,l,m)))
      }
    }
  }

  for (k in 1:nj) {
    for (l in 1:nj) {
      if(pkl(q2,q2,k,l,m) > 0) { 
        Sjj = Sjj + pkl(q2,q2,k,l,m)*log2(pkl(q2,q2,k,l,m)/(pk(q2,k,m)*pk(q2,l,m)))
      }
    }
  }

  list(NormMI=(2*Sij)/(Sii+Sjj), Sij=Sij, Sii=Sii, Sjj=Sjj)
}


matrixNormMI <- function(l) {
  n <- length(l);
  m1 <- matrix(data=0, ncol=n, nrow=n)
  m2 <- matrix(data=0, ncol=n, nrow=n)
  for (i in 1:n) {
    for (j in 1:n) {
      m1[i,j] <- exprmi(l[[i]],l[[j]])[[1]]
      m2[i,j] <- mi(l[[i]],l[[j]])[[1]]
    }
  }
  list(m1,m2)
}


