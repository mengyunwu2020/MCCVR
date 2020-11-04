library("PMA")
library("glmnet")
library(Rcpp)
library(nnet)
library("MASS")

sourceCpp(file="C:/.../MCCVRsolver.cpp")    #file = the address of MCCVRsolver.cpp


#Y is four-class response variable
#Input W,Xlist,Y to refit alpha and beta, then use Xtestlist to predict the probability on the test set.
phat<-function(W,Xlist,Y,Xtestlist)
{
  Y<-as.matrix(Y)
  K<-length(W)
  n<-nrow(Xtestlist[[1]])
  px1_test<-matrix(nrow=n,ncol=K)
  px2_test<-matrix(nrow=n,ncol=K)
  px3_test<-matrix(nrow=n,ncol=K)
  px4_test<-matrix(nrow=n,ncol=K)
  
  for (i in 1:K)
  {
    XW<-Xlist[[i]]%*%W[[i]]
    fit<-multinom(Y~XW,trace=FALSE)
    alpha1<-0
    beta1<-rep(0,ncol(W[[i]]))
    alpha2<-coef(fit)[1,1]
    beta2<-coef(fit)[1,-1]
    alpha3<-coef(fit)[2,1]
    beta3<-coef(fit)[2,-1]
    alpha4<-coef(fit)[3,1]
    beta4<-coef(fit)[3,-1]

    px1_test[,i]<-1/(1+exp(alpha2+Xtestlist[[i]]%*%W[[i]]%*%beta2-alpha1-Xtestlist[[i]]%*%W[[i]]%*%beta1)+
                       exp(alpha3+Xtestlist[[i]]%*%W[[i]]%*%beta3-alpha1-Xtestlist[[i]]%*%W[[i]]%*%beta1)+
                       exp(alpha4+Xtestlist[[i]]%*%W[[i]]%*%beta4-alpha1-Xtestlist[[i]]%*%W[[i]]%*%beta1))
    px2_test[,i]<-1/(1+exp(alpha1+Xtestlist[[i]]%*%W[[i]]%*%beta1-alpha2-Xtestlist[[i]]%*%W[[i]]%*%beta2)+
                       exp(alpha3+Xtestlist[[i]]%*%W[[i]]%*%beta3-alpha2-Xtestlist[[i]]%*%W[[i]]%*%beta2)+
                       exp(alpha4+Xtestlist[[i]]%*%W[[i]]%*%beta4-alpha2-Xtestlist[[i]]%*%W[[i]]%*%beta2))
    px3_test[,i]<-1/(1+exp(alpha1+Xtestlist[[i]]%*%W[[i]]%*%beta1-alpha3-Xtestlist[[i]]%*%W[[i]]%*%beta3)+
                       exp(alpha2+Xtestlist[[i]]%*%W[[i]]%*%beta2-alpha3-Xtestlist[[i]]%*%W[[i]]%*%beta3)+
                       exp(alpha4+Xtestlist[[i]]%*%W[[i]]%*%beta4-alpha3-Xtestlist[[i]]%*%W[[i]]%*%beta3))
    px4_test[,i]<-1/(1+exp(alpha1+Xtestlist[[i]]%*%W[[i]]%*%beta1-alpha4-Xtestlist[[i]]%*%W[[i]]%*%beta4)+
                       exp(alpha2+Xtestlist[[i]]%*%W[[i]]%*%beta2-alpha4-Xtestlist[[i]]%*%W[[i]]%*%beta4)+
                       exp(alpha3+Xtestlist[[i]]%*%W[[i]]%*%beta3-alpha4-Xtestlist[[i]]%*%W[[i]]%*%beta4))
  }
  
  px1<-apply(px1_test,1,mean)
  px2<-apply(px2_test,1,mean)
  px3<-apply(px3_test,1,mean)
  px4<-apply(px4_test,1,mean)
  
  return(list(px1=px1,px2=px2,px3=px3,px4=px4))
}




###Compute AUC in the binomial case
CalcROC <- function(y, px, ifplot = 0){ 
  ## predict binomial response with given fitted prob px, 
  ## input: 
  ## y: real response
  ## px: fitted prob seq

  n <- length(px);
  np <- sum(y == 1);
  nn <- n - np;
  sens <- c();
  spec <- c();
  
  if (np == 0) {
    auc <- sum(px < 0.5) / nn;
  } else if (nn == 0) {
    auc <- sum(px > 0.5) / np;
  } else { 
    for (i in 1:n) {
      sens[i] <- (sum(px > px[i] & y == 1) + 0.5 * sum(px == px[i] & y == 1)) / np;
      spec[i] <- (sum(px < px[i] & y == 0) + 0.5 * sum(px == px[i] & y == 0)) / nn;
    }
  }
  if (ifplot == 1) plot(1 - spec, sens);

  xx <- sort(1 - spec, index.return = TRUE)
  yy <- c(0, sens[xx$ix], 1)
  xx <- c(0, xx$x, 1)
  auc <- sum((xx[-1] - xx[-(n + 2)]) * (yy[-1] + yy[-(n + 2)]))  / 2
  
  dev <- -2 * sum(y * log(px + 0.00001) + (1 - y) * log(1 - px + 0.0001))
  return(list(spec = spec, sens = sens, auc = auc, dev = dev))
}



#Input the fitted probability px1,px2,px3,px4, compute AUC in the multinomial case
AUC_method2<-function(Y,px1,px2,px3,px4,ifplot=0)
{
  n<-nrow(Y)
  ##Y_l:label matrix
  Y_l<-matrix(0,nrow=n,ncol=4)
  for (i in 1:n)
  {
    if (Y[i,1]==1)
      Y_l[i,1]=1
    else if(Y[i,1]==2)
      Y_l[i,2]=1
    else if(Y[i,1]==3)
      Y_l[i,3]=1
    else
      Y_l[i,4]=1
  }
  ##expand the label matrix by row
  Y_lr<-matrix(nrow=1,ncol=4*n)
  for(i in 1:n)
  {
    Y_lr[1,4*i-3]<-Y_l[i,1]
    Y_lr[1,4*i-2]<-Y_l[i,2]
    Y_lr[1,4*i-1]<-Y_l[i,3]
    Y_lr[1,4*i]<-Y_l[i,4]
  }
  Y_lr<-t(Y_lr)
  
  ##px:prob matrix
  px<-cbind(px1,px2,px3,px4)
  ##expand the prob matrix by row
  px_r<-matrix(nrow=1,ncol=4*n)
  for(i in 1:n)
  {
    px_r[1,4*i-3]<-px[i,1]
    px_r[1,4*i-2]<-px[i,2]
    px_r[1,4*i-1]<-px[i,3]
    px_r[1,4*i]<-px[i,4]
  }
  px_r<-t(px_r)
  
  result<-CalcROC(Y_lr, px_r, ifplot = ifplot)
  auc<-result$auc
  return(auc)
}



#Input the fitted probability and the true Y value, compute accuracy
accuracy<-function(Y,px1,px2,px3,px4)
{
  n<-length(Y)
  Y<-as.matrix(Y)
  
  phat<-cbind(px1,px2,px3,px4)
  
  pred<-apply(phat,1,which.max)
  pred<-as.matrix(pred)
  
  
  true<-0
  for (i in 1:n)
  {
    if (Y[i,1]==pred[i,1])
      true<-true+1
  }
  accu<-true/n
  return(accu)
}

#Input the fitted probability and the true Y value, compute deviance
deviance<-function(Y,px1,px2,px3,px4)
{
  n<-length(Y)
  y<-matrix(0,nrow=n,ncol=4)
  for (i in 1:n)
    for (j in 1:4)
    {
      if (Y[i]==1)
        y[i,1]<-1
      if (Y[i]==2)
        y[i,2]<-1
      if (Y[i]==3)
        y[i,3]<-1
      if (Y[i]==4)
        y[i,4]<-1
    }
  phat<-cbind(px1,px2,px3,px4)
  
  
  sum<-0
  for (i in 1:n)
    for (j in 1:4)
    {
      add<-(-2)*y[i,j]*log(phat[i,j]+0.00001)
      sum<-sum+add
    }
  sum<-sum/n
  return(sum)
}


#generate four-class simulation data
SimulateCVR3 <- function(family = c("gaussian", "binomial", "poisson","multinomial"), n = 100, rank = 4, 
                         p1 = 50, p2 = 70, pnz = 10, sigmax = 0.2, sigmay = 0.5,   
                         beta1 = c(3, 1, 0, 0), beta2=c(1,4,0,0), beta3=c(-4,-2,0,0), beta4=c(-1,-3,0,0), standardization = TRUE){
  family <- match.arg(family)
  U <- matrix(rnorm(2 * n * rank), 2 * n, rank);
  
  V1 <- matrix(0, p1, rank);
  V11 <- matrix(runif(rank * pnz) * (rbinom(rank * pnz, 1, 0.5) * 2-1), ncol = rank);
  V2 <- matrix(0, p2, rank);
  V21 <- matrix(runif(rank * pnz) * (rbinom(rank * pnz, 1, 0.5) * 2 - 1), ncol = rank);
  V1[1:pnz, ] <- svd(V11)$u;
  V2[1:pnz, ] <- svd(V21)$u; 
  
  x1 <- U %*% t(V1); 
  x2 <- U %*% t(V2);
  lp1 <- U %*% beta1;
  lp2 <- U %*% beta2;
  lp3 <- U %*% beta3;
  lp4 <- U %*% beta4;
  e1 <- matrix(rnorm(2 * n * p1), 2 * n, p1) * sigmax;
  e2 <- matrix(rnorm(2 * n * p2), 2 * n, p2) * sigmax;
  ey1 <- rnorm(2 * n) * sigmay
  ey2 <- rnorm(2 * n) * sigmay
  ey3 <- rnorm(2 * n) * sigmay
  ey4 <- rnorm(2 * n) * sigmay
  
  lp1 <- lp1 + ey1;
  lp2 <- lp2 + ey2;
  lp3 <- lp3 + ey3;
  lp4 <- lp4 + ey4;
  
  X1 <- x1[1:n, ] + e1[1:n, ];
  X2 <- x2[1:n, ] + e2[1:n, ];
  X1test <- x1[(n + 1):(2 * n), ] + e1[(n + 1):(2 * n), ];
  X2test <- x2[(n + 1):(2 * n), ] + e2[(n + 1):(2 * n), ];
  U <- U[1:n, ]
  
  if (family == "multinomial")
  {
    px1 <- exp(lp1) / (exp(lp1)+exp(lp2)+exp(lp3)+exp(lp4));
    px2 <- exp(lp2) / (exp(lp1)+exp(lp2)+exp(lp3)+exp(lp4));
    px3 <- exp(lp3) / (exp(lp1)+exp(lp2)+exp(lp3)+exp(lp4));
    px4 <- exp(lp4) / (exp(lp1)+exp(lp2)+exp(lp3)+exp(lp4));
    y<-matrix(nrow=n,ncol=1)
    ytest<-matrix(nrow=n,ncol=1)
    for (i in c(1:n))
    {
      y[i,1]<-sample(1:4, 1, replace=T, prob = c(px1[i],px2[i],px3[i],px4[i]));
      ytest[i]<-sample(1:4, 1, replace=T, prob = c(px1[i+n],px2[i+n],px3[i+n],px4[i+n]));
    }
    if (standardization) 
    {
      V1 <- diag(apply(X1, 2, sd)) %*% V1 / sqrt(n - 1); 
      V2 <- diag(apply(X2, 2, sd)) %*% V2 / sqrt(n - 1); 
      X1 <- scale(X1); 
      X2 <- scale(X2);
      X1test <- scale(X1test); 
      X2test <- scale(X2test);
    }
    snr.x1 <- sd(x1) / sd(e1)
    snr.x2 <- sd(x2) / sd(e2)
    snr <- c(snr.x1, snr.x2, 0)
    
  }
  return(list(X1 = X1, X2 = X2, y = as.matrix(y), U = U, V1 = V1, V2 = V2,X1test=X1test, X2test=X2test,ytest=ytest))
}




#use CCA to initialize W, when there are two X datasets.
SparseCCA <- function(X1, X2, rank = 2, penaltyx1 = NULL, penaltyx2 = NULL, nperms = 25, ifplot = 0){
  if (nrow(X1) != nrow(X2)) stop("X1 and X2 should have the same row number!")
  if (ncol(X1) == 1 | ncol(X2) == 1) stop("X1 and X2 should have at least 2 columns!")
  if (is.null(penaltyx1) | is.null(penaltyx2)){
    penaltyx1 <- seq(0.1, 0.7, len = 20)
    penaltyx2 <- seq(0.1, 0.7, len = 20)
  }
  if (length(penaltyx1) > 1){
    perm.out <- CCA.permute(X1, X2, "standard", "standard", penaltyx1, penaltyx2, nperms = nperms)  
    if (ifplot == 1 ) {
      print(c(perm.out$bestpenaltyx, perm.out$bestpenaltyz))
      plot(perm.out)
    }
    SparseCCA_X12 <- CCA(X1, X2, typex = "standard", typez = "standard", K = rank,
                         penaltyx = perm.out$bestpenaltyx, penaltyz = perm.out$bestpenaltyz,
                         v = perm.out$v.init, trace = FALSE)
  } else {
    SparseCCA_X12 <- CCA(X1, X2, typex = "standard", typez = "standard", K = rank, 
                         penaltyx = penaltyx1, penaltyz = penaltyx2, trace = FALSE)
  }
  D1 <- diag(diag(t(X1 %*% SparseCCA_X12$u) %*% X1 %*% SparseCCA_X12$u)^(-0.5), rank, rank)
  D2 <- diag(diag(t(X2 %*% SparseCCA_X12$v) %*% X2 %*% SparseCCA_X12$v)^(-0.5), rank, rank)
  W1 <- SparseCCA_X12$u %*% D1
  W2 <- SparseCCA_X12$v %*% D2 ## to make X1 %*% W1 orthogonal
  return(list(W1 = W1, W2 = W2))
}


#use CCA to initialize W, when there are more than two X datasets.
SparseCCA_xlist <- function(Xlist, rank = 2, penalties=NULL){
  K<-length(Xlist)
  n<-nrow(Xlist[[1]])
  for (i in 1:K)
  {
    if (nrow(Xlist[[i]])!=n)
      stop("All Xs should have the same row number!")
  }

  for (i in 1:K)
  {
    if (ncol(Xlist[[i]])==1)
      stop("All Xs should have at least 2 columns!")
  }
  minpk<-ncol(Xlist[[1]])
  for (i in 1:K)
  {
    if (ncol(Xlist[[i]])<minpk)
      minpk<-ncol(Xlist[[i]])
  }
 
  if (is.null(penalties)){
    penalties <- seq(1.1, sqrt(minpk), len = 20)
  }
  if (length(penalties) > 1)
  {
    multiperm.out <- MultiCCA.permute(Xlist, type="standard",penalties = penalties,trace=FALSE)
    multiSparseCCA_X12 <- MultiCCA(Xlist, penalty=multiperm.out$bestpenalties, ws=multiperm.out$ws.init, type="standard",  ncomponents=rank, trace=FALSE)
  } else 
  {
    multiSparseCCA_X12<-MultiCCA(Xlist,penalty=penalties, type="standard",ncomponents = rank,trace=FALSE)
  }
  D<-list()
  W<-list()
  for (i in 1:K)
  {
    D[[i]]<-diag(diag(t(Xlist[[i]] %*% multiSparseCCA_X12$ws[[i]]) %*% Xlist[[i]] %*% multiSparseCCA_X12$ws[[i]])^(-0.5), rank, rank)
    W[[i]]<-multiSparseCCA_X12$ws[[i]] %*% D[[i]]
  }  ## to make X1 %*% W1 orthogonal
  
  return(W)
}




#Use cross validation to select parameters of MCCVR, when there are two X datasets.
#type.measure="auc", "accuracy", "deviance" 
multiCVR_X1X2 <- function(Y, Xlist, rankseq = 2, neta = 10, etaseq = NULL, nlam = 50,  
                          Lamseq = NULL, Wini = NULL, penalty = c("GL1", "L1"), nfold = 10,
                          foldid = NULL, opts = list(), type.measure = NULL)
{
  family <- "multinomial"
  penalty <- match.arg(penalty)
  
  if (is.null(type.measure))
    type.measure<-"auc"
  
  if(is.null(penalty)) {
    penalty = "GL1"; 
    cat('Use default: penalty = "GL1".\n')
  }
  
  Y <- as.matrix(Y)
  X1 <- as.matrix(Xlist[[1]])
  X2 <- as.matrix(Xlist[[2]])
  Xlist <- list(X1 = X1, X2 = X2)
  n <- dim(Y)[1]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  
  if(p1 == 1 | p2 == 1) 
    stop("X1 and X2 should have at least 2 columns!")
  if(nrow(X1) != nrow(X2)) 
    stop("X1 and X2 should have the same row number!")
  if(n != nrow(X1)) 
    stop("Y and X1, X2 should have the same row number!")
  
  opts$n  <- n
  opts$p1 <- p1
  opts$p2 <- p2
  
  
  
  if (is.null(opts$standardization)) {
    opts$standardization <- TRUE; 
    #cat('Use default: standardization = TRUE.\n')
  }  
  if (is.null(opts$maxIters)) {
    opts$maxIters <- 500; 
    #cat('Use default: maxIters = 500.\n')
  }
  if (is.null(opts$tol)) {
    opts$tol <- 0.01; 
    #cat('Use default: tol = 0.01.\n')
  }
  if (is.null(opts$spthresh)) {
    spthresh <- 0.4
    opts$spthresh = spthresh; 
    #cat('Use default: spthresh = 0.4.\n')
  }else{
    spthresh<-opts$spthresh
  }
    
  
  if (is.null(etaseq)) {
    if(is.null(neta)) 
      neta <- 10
    etaseq <- seq(0.01, 0.99, len = neta)
  }
  if (is.null(Lamseq)) {
    if (is.null(nlam)) 
      nlam <- 50
    lamseq <- seq(0.1, 20, len = nlam)
    Lamseq <- cbind(lamseq, lamseq*p2/p1)
  }
  if (is.null(foldid)) { 
    if(is.null(nfold))
      nfold <- 10
    tmp <- ceiling(c(1:n)/(n/nfold));
    foldid <- sample(tmp, n) 
  } else {
    if(nfold != length(unique(foldid))){
      nfold <- length(unique(foldid)) 
    }
  }
  
  
  
  Lr <- length(rankseq)
  Leta <- length(etaseq)
  Llam<-nrow(Lamseq)
  
  pred<-array(0,c(nfold,Lr,Leta,Llam))
  sparse<-array(0,c(nfold,Lr,Leta,Llam))
  
  warm<-FALSE
  
  for (f in 1:nfold)
  {
    X1train <- Xlist[[1]][foldid != f, ];
    X2train <- Xlist[[2]][foldid != f, ];
    Xtrain  <- list(X1 = X1train, X2 = X2train);
    Ytrain  <- as.matrix(Y[foldid != f,  ]);
    
    X1test  <- Xlist[[1]][foldid == f, ];
    X2test  <- Xlist[[2]][foldid == f, ];
    Xtest   <- list(X1test, X2test);
    Ytest   <- as.matrix(Y[foldid == f, ]);
    
    for (a in 1:Lr)
    {
      Wini<-SparseCCA(X1, X2, rankseq[a], 0.7, 0.7);
      Wwork<-Wini
      for (b in 1:Leta)
        for (c in 1:Llam)
        {
          r<-rankseq[a]
          eta<-etaseq[b]
          lam<-Lamseq[c,]
          
          if (c == 1)
          {
            optswork <- opts;
            optswork$maxIters <- 2 * opts$maxIters;
          } else {
            optswork <- opts;
            if (warm)
              Wwork <- fitwork$W; ## warm start
          }
          
          fitwork <- multicvrsolver(Ytrain, Xtrain, r, eta, as.vector(lam), family, Wwork, penalty, optswork)
          W1<-fitwork$W[[1]]
          W2<-fitwork$W[[2]]
          
          p<-phat(fitwork$W, Xtrain, Ytrain, Xtest)
          if (type.measure=="auc"){
            pred[f,a,b,c]<-AUC_method2(Ytest, p$px1, p$px2, p$px3, p$px4, ifplot=0)
          }else if(type.measure=="accuracy"){
            pred[f,a,b,c]<-accuracy(Ytest, p$px1, p$px2, p$px3, p$px4)
          }else if(type.measure=="deviance"){
            pred[f,a,b,c]<-(-1)*deviance(Ytest, p$px1, p$px2, p$px3, p$px4)
          }
          
          sparse[f,a,b,c]<-(sum(W1 != 0) + sum(W2 != 0)) / (length(W1) + length(W2))
        }  #loop of Llam
    }  #loop of Lr
  }  #loop of nfold
  
  maxpred<- (-100000)
  ir<-0
  ieta<-0
  ilam<-0
  for (a in 1:Lr)
    for (b in 1:Leta)
      for (c in 1:Llam)
      {
        if (mean(pred[,a,b,c])>maxpred & mean(sparse[,a,b,c])<spthresh)
        {
          maxpred<-mean(pred[,a,b,c])
          ir<-a
          ieta<-b
          ilam<-c
        }
      }
  rankhat<-rankseq[ir]
  etahat<-etaseq[ieta]
  lamhat<-Lamseq[ilam,]
  
  Winihat <- SparseCCA(X1, X2, rankhat, 0.7, 0.7)
  fit<-multicvrsolver(Y, Xlist, rankhat, etahat, as.vector(lamhat), family, Winihat, penalty, opts)
  
  return(list(fit=fit,rankhat=rankhat,etahat=etahat,lamhat=lamhat,pred=pred,sparse=sparse,maxpred=maxpred,pred=pred,ir=ir,ieta=ieta,ilam=ilam))
}


#Use cross validation to select parameters of MCCVR, when there are more than two X datasets.
multiCVR_Xlist <- function(Y, Xlist, rankseq = 2, neta = 10, etaseq = NULL, nlam = 50,  
                           Lamseq = NULL, Wini = NULL, penalty = c("GL1", "L1"), nfold = 10,
                           foldid = NULL, opts = list(), type.measure = NULL)
{
  family <- "multinomial"
  penalty <- match.arg(penalty)
  
  if (is.null(type.measure))
    type.measure<-"auc"
  
  if(is.null(penalty)) {
    penalty = "GL1"; 
    cat('Use default: penalty = "GL1".\n')
  }
  
  K<-length(Xlist)
  
  Y <- as.matrix(Y)
  n <- dim(Y)[1]
  p<-c()
  for (k in 1:K)
  {
    p[k] <- dim(Xlist[[k]])[2]
    if(p[k]==1)
      stop("All Xs should have at least 2 columns!")
    if(nrow(Xlist[[k]]) != n)
      stop("Y and All Xs should have the same row number!")
  }
  
  opts$n  <- n
  opts$p  <- p
  
  
  minpk<-ncol(Xlist[[1]])
  for (i in 1:K)
  {
    if (ncol(Xlist[[i]])<minpk)
      minpk<-ncol(Xlist[[i]])
  }
  
  
  if (is.null(opts$standardization)) {
    opts$standardization <- TRUE; 
    #cat('Use default: standardization = TRUE.\n')
  }  
  if (is.null(opts$maxIters)) {
    opts$maxIters <- 500; 
    #cat('Use default: maxIters = 500.\n')
  }
  if (is.null(opts$tol)) {
    opts$tol <- 0.01; 
    #cat('Use default: tol = 0.01.\n')
  }
  if (is.null(opts$spthresh)) {
    spthresh <- 0.4
    opts$spthresh = spthresh; 
    #cat('Use default: spthresh = 0.4.\n')
  }else{
    spthresh<-opts$spthresh
  } 
  
  
  
  if (is.null(etaseq)){
    if(is.null(neta)) 
      neta <- 10
    etaseq <- seq(0.01, 0.99, len = neta)
  }
  if (is.null(Lamseq)){
    if (is.null(nlam)) 
      nlam <- 50
    lamseq <- seq(0.1, 20, len = nlam)
    Lamseq<-lamseq
    for (j in 2:K)
    {
      Lamseq<-cbind(Lamseq,lamseq*p[j]/p[1])
    }
  } 
  
  if (is.null(foldid)){ 
    if(is.null(nfold))
      nfold <- 10
    tmp <- ceiling(c(1:n)/(n/nfold));
    foldid <- sample(tmp, n) 
  }else{
    if(nfold != length(unique(foldid))){
      nfold <- length(unique(foldid)) 
    }
  }

  Lr <- length(rankseq)
  Leta <- length(etaseq)
  Llam<-nrow(Lamseq)
  
  pred<-array(0,c(nfold,Lr,Leta,Llam))
  sparse<-array(0,c(nfold,Lr,Leta,Llam))
  
  warm<-FALSE
  
  for (f in 1:nfold)
  {
    Xtrain<-list()
    Xtest<-list()
    Ytrain  <- as.matrix(Y[foldid != f,  ])
    Ytest   <- as.matrix(Y[foldid == f, ])
    for (t in 1:K)
    {
      Xtrain[[t]]<-Xlist[[t]][foldid != f,]
      Xtest[[t]]<-Xlist[[t]][foldid==f,]
    }
    
    for (a in 1:Lr)
    {
      Wini<-SparseCCA_xlist(Xlist, rankseq[a], max(0.7*sqrt(minpk),1.1))
      Wwork<-Wini
      for (b in 1:Leta)
        for (c in 1:Llam)
        {
          r<-rankseq[a]
          eta<-etaseq[b]
          lam<-Lamseq[c,]
          
          if (c == 1)
          {
            optswork <- opts;
            optswork$maxIters <- 2 * opts$maxIters;
          } else {
            optswork <- opts;
            if (warm)
              Wwork <- fitwork$W; ## warm start
          }
          
          
          fitwork <- multicvrsolver(Ytrain, Xtrain, r, eta, as.vector(lam), family, Wwork, penalty, optswork)
          Wlist<-fitwork$W
          
          p<-phat(fitwork$W, Xtrain, Ytrain, Xtest)
          if (type.measure=="auc"){
            pred[f,a,b,c]<-AUC_method2(Ytest, p$px1, p$px2, p$px3, p$px4, ifplot=0)
          }else if(type.measure=="accuracy"){
            pred[f,a,b,c]<-accuracy(Ytest, p$px1, p$px2, p$px3, p$px4)
          }else if(type.measure=="deviance"){
            pred[f,a,b,c]<-(-1)*deviance(Ytest, p$px1, p$px2, p$px3, p$px4)
          }
          
          sum_nonzero<-0
          sum_length<-0
          for (v in 1:K)
          {
            sum_nonzero<-sum_nonzero+sum(Wlist[[v]]!=0)
            sum_length<-sum_length+sum(length(Wlist[[v]]))
          }
          ## sparse = proportion of nonzero elements among W1 and W2
          sparse[f,a,b,c] <- sum_nonzero / sum_length  
        }
    }
    
  }
  
  
  
  maxpred<- (-100000)
  ir<-0
  ieta<-0
  ilam<-0
  for (a in 1:Lr)
    for (b in 1:Leta)
      for (c in 1:Llam)
      {
        if (mean(pred[,a,b,c])>maxpred & mean(sparse[,a,b,c])<spthresh)
        {
          maxpred<-mean(pred[,a,b,c])
          ir<-a
          ieta<-b
          ilam<-c
        }
      }
  rankhat<-rankseq[ir]
  etahat<-etaseq[ieta]
  lamhat<-Lamseq[ilam,]
  
  Winihat <- SparseCCA_xlist(Xlist, rankhat, max(0.7*sqrt(minpk),1.1))
  fit<-multicvrsolver(Y, Xlist, rankhat, etahat, as.vector(lamhat), family, Winihat, penalty, opts)
  
  return(list(fit=fit,rankhat=rankhat,etahat=etahat,lamhat=lamhat))
}




