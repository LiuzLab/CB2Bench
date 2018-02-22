fit.ab2 <- function(nvec,xvec){
  # FIT.AB2(NVEC,XVEC) -- Given the vectors of library sizes
  # (nvec) and associated counts for a particular tag (xvec),
  # this fits the alpha and beta values of the underlying beta
  # distribution using the method of moments. It returns
  # PHAT - the estimated proportion
  # VHAT - the associated variance estimate
  # VMIN - the minimum possible variance
  # ALPHAHAT - the first of the beta distribution parameters
  # BETAHAT - the second of the beta distribution parameters

  # This script was written to accompany
  # "Differential Expression in SAGE: Accounting for
  # Normal Between-Library Variation" by Baggerly,
  # Deng, Morris and Aldaz, Bioinformatics, 2003.
  #
  # Keith Baggerly, kabagg@mdanderson.org
  # Last Updated Dec 5, 2002

  n <- length(nvec);

  # Get the associated proportions

  pvec <- xvec / nvec;

  # assume a starting weight vector proportional to
  # the library sizes

  wvec <- nvec / sum(nvec);

  # Get the minimum variance

  phat <- sum(wvec * pvec);
  vmin <- phat*(1-phat)/sum(nvec);

  # Iterate to convergence (this is generally quite rapid)

  nrep <- 15;

  resmat1 <- matrix(0,nrep,4);
  resmat2 <- matrix(0,nrep,n);
  kvec <- rep(0,2);
  kvec2 <- rep(0,4);
  counter <- 1;

  for(i in 1:nrep){

    phat <- sum(wvec * pvec);
    w2 <- sum(wvec * wvec);
    vhat <- (sum(wvec^2 * pvec^2) - w2*(phat^2)) / (1 - w2);

    if(vhat > vmin){

      # Fit beta and alpha

      w2n <- sum(wvec^2 / nvec);
      betahat <- (vhat - phat*(1-phat)*w2) / (phat*w2n - (vhat/(1-phat)));
      alphahat <- (phat/(1-phat))*betahat;

      # refit the weight vector

      wvec <- ((alphahat + betahat)*nvec) / (alphahat + betahat + nvec);
      wvec <- wvec / sum(wvec);

    }else{

      alphahat <- 0;
      betahat <- 0;
      vhat <- vmin;
      break;

    }

    resmat1[i,] <- c(phat, vhat, betahat, alphahat);
    resmat2[i,] <- wvec;

    # Tack on a bisection approach in case convergence
    # seems to be slow

    kvec[((i+1) %% 2) + 1] = alphahat + betahat;
    if(!(i %% 2)){
      ktemp <- mean(kvec);
      wvec <- (ktemp * nvec) / (ktemp + nvec);
      wvec <- wvec / sum(wvec);
    }

  }

  return(list(phat=phat,vhat=vhat,vmin=vmin,alphahat=alphahat,betahat=betahat))

}

###### Baggerly t-test function #######
bbt.test = function(Amat, Bmat){
  nvecA = apply(Amat,2,sum)
  nvecB = apply(Bmat,2,sum)
  ntarget = nrow(Amat)
  nA = ncol(Amat)
  nB = ncol(Bmat)
  stat = matrix(NA,nrow=ntarget,ncol=9)
  colnames(stat) = c("tstat","df","pvalue.twoside","pvalue.pos", "pvalue.neg", "propA","varA","propB","varB")
  for(i in 1:ntarget){
    fitA = fit.ab2(nvecA,Amat[i,])
    fitB = fit.ab2(nvecB,Bmat[i,])
    zeroVar = (fitA$vhat == 0) & (fitB$vhat == 0)
    #t-statistic
    stat[i,1] = (fitA$phat - fitB$phat)/sqrt(fitA$vhat + fitB$vhat + zeroVar)
    #variance
    stat[i,2] = ((fitA$vhat + fitB$vhat)^2+zeroVar)/((fitA$vhat)^2/(nA-1) + (fitB$vhat)^2/(nB-1)+zeroVar)
    #p-value
    stat[i,3] = 2*pt(-abs(stat[i,1]),stat[i,2])
    stat[i,4] = 1-pt(stat[i,1],stat[i,2])
    stat[i,5] = 1-pt(stat[i,1],stat[i,2], lower.tail = F)

    stat[i,6] = fitA$phat
    stat[i,7] = fitA$vhat
    stat[i,8] = fitB$phat
    stat[i,9] = fitB$vhat
  }
  stat
}

