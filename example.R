sim_beta1 = c(3,2,0,0)
sim_beta2=c(1,4,0,0)
sim_beta3=c(-4,-1,0,0)
sim_beta4=c(-2,-3,0,0)

#generate simulation data
mydata<-SimulateCVR3(family="multinomial",n=100,sigmax = 0.2, sigmay = 0.5,p1=50,p2=70,beta1=sim_beta1,beta2=sim_beta2,beta3=sim_beta3,beta4=sim_beta4)
X1 <- mydata$X1
X2 <- mydata$X2
Xlist <- list(X1 = X1, X2 = X2)
Y <- mydata$y

X1test<-mydata$X1test
X2test<-mydata$X2test
Xtestlist<-list(X1test,X2test)
Ytest<-mydata$ytest


Lamseq<-10^seq(-1,0.7,len=10)
Lamseq<-cbind(Lamseq,Lamseq)
etaseq<-0.1*c(1:9)
opts<-list(standardization = TRUE, maxIters = 300, tol = 0.01,spthresh=0.99)
#fit the model
fit<-multiCVR_X1X2(Y, Xlist, rankseq = 4, neta = length(etaseq), etaseq = etaseq, nlam = nrow(Lamseq),
                   Lamseq = Lamseq, Wini = NULL, penalty = "GL1", nfold = 5,
                   foldid = NULL, opts = opts, type.measure="accuracy")
W1<-fit$fit$W[[1]]
W2<-fit$fit$W[[2]]
W <- list(W1,W2)

prob <- phat(W,Xlist,Y,Xtestlist)
AUC_method2(Ytest,prob$px1,prob$px2,prob$px3,prob$px4,ifplot=1)    #compute AUC
accuracy(Ytest,prob$px1,prob$px2,prob$px3,prob$px4)    #compute accuracy
