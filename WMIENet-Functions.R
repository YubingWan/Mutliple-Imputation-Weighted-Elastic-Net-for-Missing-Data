#############################################################################
library(mice)
#############################################################################
##########        Self-defined Preprocessing Funtions          ##############
#############################################################################
## Function to scale and standardize covariate matrix X
my.stdz <- function (x, center = TRUE, scale = TRUE, ... ) 
{
    nc <- ncol(x)
    if (center) {
        means <- apply(x, 2, mean, na.rm=TRUE)
        x     <- sweep(x, 2, means)
    }
    if (scale) {
        sd <- apply(x, 2, sd, na.rm=TRUE)
        x  <- sweep(x, 2, sd, "/")
    }
    x 
}

## Function to create weight vector for each feature of covariate matrix X
## M is the total number of multiple imputation to perform
wtVar <- function(x, M) 
{
   m   <- dim(x)[2]  
   wts <- NULL
   for (i in 1:m){ 
       f   <- length(x[is.na(x[,i]), i])
       s   <- length(x[,i])
       wts <- c(wts, (1-f/s)/M)} 
   return(wts)
}

## Another version of function to create weight vector for each feature of covariate matrix X
## M is the total number of multiple imputation to perform
wtObs <- function(x, M) 
{
   m   <- dim(x)[1]  
   wts <- NULL
   for (i in 1:m){ 
       f   <- length(x[i, is.na(x[i,])])
       s   <- length(x[i,])
       wts <- c(wts, (1-f/s)/M)} 
   wts <- rep(wts, M)
   return(wts)
}

####################################################################################
##########  Functions for Estimating Penalized Regression Coefficients  ############
####################################################################################
## Thresholding function 
S.wnet <- function(z, r, ...)
{  result<-0
   if (z>0 & r<abs(z)) result <- z-r
   if (z<0 & r<abs(z)) result <- z+r
   return(result)
}

## Function to solve for regression coefficient estimate of weighted elastic net model 
## for the multiple-imputed data
Solve.Beta.WENET <- function(X, Y, W, M=1, lambda, alpha, ...)
{
   N      <- nrow(X)/M
   if (missing(W)) {W = rep(1/M, N)}
   X1.new <- as.matrix(cbind(rep(1, nrow(X)), X))
   XWY    <- t(X1.new)%*%(W*Y)/N
   XWX    <- t(X1.new)%*%diag(W)%*%X1.new/N
   beta0  <- beta.new <- rep(0, length(X1.new[1, ]))
   beta.new[length(beta0)] <- beta0[length(beta0)]+1
   while(max(abs(beta.new-beta0)) > 0.0001)
   {  beta0 <- beta.new
      for (j in 1:length(beta0))
      {  z1 <- XWY[j]
         if(j > 1) {z1 <- z1-XWX[j,1:(j-1)]%*%beta.new[1:(j-1)]}
         if(j < length(beta0)) {z1 <- z1-XWX[j,(j+1):length(beta0)]%*%beta0[(j+1):length(beta0)]}
         beta.new[j] <- S.wnet(z1, lambda*alpha)/(XWX[j,j]+lambda*(1-alpha))
      }
   }
   beta.new.undo <- as.vector((diag(XWX)+(1-alpha)*lambda)/diag(XWX))*beta.new
   return(list(beta=beta.new, beta.undo=beta.new.undo))
}

## Function to perform cross-validation to decide best estimate of beta
calculate.cv <- function(X, Y, W, M=1, lambda, alpha, nfolds=10, foldid, ...)
{  
   N <- nrow(X)/M
   if (missing(W)) {W = rep(1/M, N)}
   cv.err.undo <- 0; beta.est.undo <- c()
   if (missing(foldid)) {foldid = rep(sample(rep(seq(nfolds), length=N)), M)}
   for (i in unique(foldid)) {
       which = foldid == i
       y.sub = Y[!which]
       Result <- Solve.Beta.WENET(X=X[!which,], Y=y.sub, W=W[!which], M=M, lambda=lambda, alpha=alpha, ...)
       beta.est.undo <- rbind(beta.est.undo, Result$beta.undo)
       cv.err.undo <- cv.err.undo+sum(W[which]*(Y[which]-Result$beta.undo[1]-X[which,]%*%Result$beta.undo[c(-1)])^2)
   }
   cv.err.undo <- cv.err.undo/N
   return(list(cv.err.undo=cv.err.undo, beta.undo=beta.est.undo))
}

## Tunning function to estimate best penalized estimate of beta with function defined above
tuning.wnet <- function (X, Y, W, M=1, alphas, nfolds=10, foldid, K=100, eps=0.01, ...) 
{
    BETA.undo <- PMSE.undo <- LAMBDA <- NULL
    N   <- nrow(X)/M
    if (missing(W)) {W = rep(1/M, N)}
    X1  <- as.matrix(cbind(rep(1, nrow(X)), X))
    Y1  <- rep(Y, 1)
    XWY <- t(X1)%*%(W*Y1)/N
    if (missing(foldid)) {foldid=rep(sample(rep(seq(nfolds), length=N)), M)}  
    for (i in 1:length(alphas)) {
       alpha   <- ifelse(alphas[i] > 0, alphas[i], 0.001) 
       lam.max <- max(abs(XWY))/alpha  
       lam.seq <- exp(seq(log(lam.max*eps), log(lam.max), (log(lam.max)-log(lam.max*eps))/K))     
       my.beta=my.beta.undo=my.PMSE=my.PMSE.undo=my.PMSE1=my.PMSE.undo1=penalty=c()
       for (lambda in lam.seq) {
           result       <- Solve.Beta.WENET(X=X,Y=Y,W=W,M=M,lambda=lambda,alpha=alpha,...)
           S.beta.undo  <- result$beta.undo           
           CV.ERR <- calculate.cv(X=X,Y=Y,W=W,M=M,lambda=lambda,alpha=alpha,foldid=foldid,...)          
             # Bcv    <- CV.ERR$beta.undo
             # Bcv.mn <- apply(Bcv, 2, mean)
             # Bcv.sd <- apply(Bcv, 2, sd)
             # Bcv.er <- qt(0.975, df=nfolds-1)*Bcv.sd
             # CI     <- cbind(Bcv.mn-Bcv.er, Bcv.mn+Bcv.er)
             # for (j in 1:length(S.beta.undo)) {if (CI[j,1]*CI[j,2]<0) S.beta.undo[j]=0}
           my.beta.undo <- rbind(my.beta.undo, S.beta.undo)
           my.PMSE.undo <- c(my.PMSE.undo, CV.ERR$cv.err.undo)           
       }      
       index  <- which(my.PMSE.undo <= min(my.PMSE.undo))
       LAMBDA <- c(LAMBDA, lam.seq[index])
       BETA.undo <- rbind(BETA.undo, my.beta.undo[index, ])
       PMSE.undo <- c(PMSE.undo, my.PMSE.undo[index])
    }    
    index    <- which(PMSE.undo <= min(PMSE.undo))
    Beta.opt <- BETA.undo[index,]
    PMSE.opt <- PMSE.undo[index]
    Alph.opt <- alphas[index]
    Lamb.opt <- LAMBDA[index]   
    return(list(Beta.opt=Beta.opt, Alph.opt=Alph.opt, Lamb.opt=Lamb.opt, PMSE.opt=PMSE.opt))
}

## Self-defined summary function
sum.wnet <- function(Obj, xTest, yTest, sensInd, specInd, ...)
{ 
   coef <- Obj$Beta.opt[-1]
   pmse <- mean((yTest-Obj$Beta.opt[1]-xTest%*%Obj$Beta.opt[-1])^2)
   seVA <- coef[sensInd]
   spVA <- coef[specInd]
   sens <- length(seVA[abs(seVA)!=0])/length(sensInd)
   spec <- length(spVA[abs(spVA)==0])/length(specInd)
   return (list(COEF=coef, PMSE=pmse, SENS=sens, SPEC=spec))
} 

