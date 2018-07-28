
EMSE=function(D,C,simplecomparison=TRUE){
  N=length(D)
  N1=sum(D)
  N0=N-N1
  
  if (simplecomparison){ #EMSE for difference in means
    variance=1/N1 + 1/N0 #variance of difference in means, assuming residual variance of 1
    w1=(1/N1)*D - rep(1/N, N)
    w0= -(1/N0)*(1-D) + rep(1/N, N)
    Ebias2= w1 %*% C %*% w1 + w0 %*% C %*% w0
    MSE=variance+Ebias2
    
  } else { #EMSE for optimal estimator
    Cbar=colMeans(C)
    priorvar= 2 * mean(Cbar) #prior variance of difference in means
    v1=Cbar[D] %*% solve(C[D,D]+diag(N1), Cbar[D]) #prior variance of posterior expectation for E[Y^1]
    v0=Cbar[!D] %*% solve(C[!D,!D]+diag(N0), Cbar[!D]) #same for Y^0
    MSE=priorvar-v1-v0
  }
  
  MSE
}


#function to generate random N vector of N1 1s and N-N1 0s
binarysample=function(N, N1){
  D=rep(0,N)
  D[sample(N,N1)]=1
  as.logical(D)
}

#function to rerandomize R times, and find optimal assignment minimizing EMSE
maxEMSE=function(C, R=1000,
                 simplecomparison=TRUE,
                 parallel=TRUE){
  #number of units and number of treated units
  N=dim(C)[1]
  N1=floor(N/2)
  
  #matrix with binary column vectors of treatment assignments
  Ds=sapply(1:R, function(i) binarysample(N,N1))
  
  #vector of corresponding expected MSEs (in parallel, if so desired)
  if (parallel) {
    library(parallel)
    no_cores = detectCores()
    clust = makeCluster(no_cores, type="FORK")  #forking requires mac or Linux OS!
    EMSEs=parSapply(clust, 1:R, function(i) EMSE(Ds[,i], C, simplecomparison=simplecomparison))
    stopCluster(clust)
  } else {
    EMSEs=sapply(1:R, function(i) EMSE(Ds[,i], C, simplecomparison=simplecomparison))
  }
  
  #assignment which minimizes EMSEs
  i_optimal_D=which.min(EMSEs)
  list(Dstar=as.integer(Ds[,i_optimal_D]),
       bestMSE=EMSEs[i_optimal_D],
       meanMSE=mean(EMSEs))
}





#prior covariance matrix, relative to an assumed residual variance of 1
Csquaredexponential=function(covariates, #dataframe of covariates
                             smootW=.2,  #smoothness of prior: smaller = smoother; inverse of lengthscale
                             R2=.5 #prior expectation of R2 of covariates, within treatment arm
                             ){

  covariates=as.data.frame(covariates)
  N=dim(covariates)[1] #number of observations
 
  #distance squared betwee covariates
  distance2=matrix(0, N, N)
  for (j in names(covariates)){
     distance2 = distance2 +(smootW^2/var(covariates[[j]]))*outer(covariates[[j]], covariates[[j]], "-")^2 #add distance squared, normalized by std and lengthscale
  }
  #squared exponential covariance kernel
  C=exp(-distance2/2)
  
  #rescale C to get required R2
  predictionvar=mean(diag(C)) - mean(C)/N
  rescale=(R2/(1-R2))/predictionvar
  rescale*C
}

Clinear=function(covariates, #dataframe of covariates
                 Sigma=diag(dim(covariates)[2]+1), #variance matrix of coefficients and intercept
                 R2=.5 #prior expectation of R2 of covariates, within treatment arm
                 ){
  N=dim(covariates)[1] #number of observations
  X=  as.matrix(covariates)
  X= X - rep(colMeans(X), rep.int(nrow(X), ncol(X))) #demean columns
  X=cbind(rep(1,N), X) #include intercept
  
  C= X %*% Sigma %*% t(X)
  
  #rescale C to get required R2
  predictionvar=mean(diag(C)) - mean(C)/N
  rescale=(R2/(1-R2))/predictionvar
  rescale*C
}




C_from_covariates_file=function(filename,
                                header=FALSE,
                                squaredexp=TRUE
                                ){
  covariates=read.table(filename,header=header, sep=",")
  
  if (squaredexp) Csquaredexponential(covariates)
    else Clinear(covariates)
  
}


