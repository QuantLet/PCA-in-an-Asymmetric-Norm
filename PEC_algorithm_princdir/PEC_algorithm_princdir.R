# ------------------------------------------------------------------------------
# Book:       	N.M. Tran, M. Osipenko, and W.K. Härdle 
#               "Principal Components in an Asymmetric Norm"
# ------------------------------------------------------------------------------
# Quantlet:   	PrincipalDirection
# ------------------------------------------------------------------------------
# Description:	Functional library for the PrincipalExpectile algorithm
# ------------------------------------------------------------------------------
# See also:     China_example.R, fpca.R
# ------------------------------------------------------------------------------
# Keywords:     PCA, principal components, expectile, quantile
# ------------------------------------------------------------------------------
# Usage:      	-
# ------------------------------------------------------------------------------
# Inputs:     	-
# ------------------------------------------------------------------------------
# Output:     	-
# ------------------------------------------------------------------------------
# Example:      -
# ------------------------------------------------------------------------------
# Author:     	Ngoc Mai Tran, Maria Osipenko
# ------------------------------------------------------------------------------

#-- Functional library of the PrincipalExpectile algorithm

#FUNCTION mhat
#   return the tau-expectile estimate with given partition
#   Y: data matrix p times n. p is the dimension, n number of observations
#   lab: label vector n times 1. Specifies the partitions. 1 for positive, -1 for negative.
#   alpha: -0.5 < alpha < 0.5. alpha + 0.5 is the expectile level tau.
#output: 

mhat <- function(Y, lab, alpha){
    p = dim(Y)[1]
    nplus <- sum(lab == 1)
    nminus <- sum(lab == -1)
    if(p > 1){
        sumplus <- apply(as.matrix(Y[,lab == 1]), 1, sum)
        summinus <- apply(as.matrix(Y[,lab == -1]), 1, sum)
    } else{
        sumplus <- sum(Y[,lab == 1])
        summinus <- sum(Y[,lab == -1])
    }
    top <- (0.5 + alpha)*sumplus + (0.5-alpha)*summinus
    bottom <- (0.5+alpha)*nplus + (0.5-alpha)*nminus
    return(top/bottom)
}

#FUNCTION expectile. iteratively find the tau-expectile of a sequence of numbers in R

expectile <- function(Y, alpha,nurexp=TRUE){
    Y=matrix(Y,1,length(Y))
    n = dim(Y)[2]
    lab <- ifelse(c(1:n) > n/2, 1, -1)        
    change = TRUE
    while(change){
        change = FALSE
        ehat <- mhat(Y, lab, alpha)
        newlab <- ifelse(Y > ehat, 1, -1)
        change <- ifelse(sum(lab!=newlab)>0, TRUE, FALSE)
        lab <- newlab
    }
	if (nurexp==FALSE){
		out=lab
	} else {
		out=ehat
	}
    return(out)
}

#DEBUG FUNCTION: compute explicitly the value of the loss function
#given Y, psi, alpha.
loss <- function(Y, psi, alpha){
    obs = t(psi)%*%Y
    ehat = expectile(obs, alpha)
    sign = ifelse(obs-ehat > 0, 0.5+alpha, 0.5-alpha)
    return(sum((obs-ehat)^2*sign))
}

#FUNCTION: return the first eigenvector psi given the partition
psihat <- function(Y, lab, alpha){
    n = dim(Y)[2]
    m <- mhat(Y, lab, alpha)
    #compute the covariance matrix normalized by m instead of the mean
    Yplus <- Y[,lab == 1]
    Yminus <- Y[,lab == -1]    
    Cplus = (0.5+alpha)/n*(Yplus-m)%*%t(Yplus-m)
    Cminus = (0.5-alpha)/n*(Yminus-m)%*%t(Yminus-m)
    C = Cplus + Cminus
    return(eigen(C)$vectors[,1])
}

#FUNCTION: given a direction psi, compute the expectile of the data in this direction. Return the lab.
expdir <- function(Y, psi, alpha){
    #compute the inner product
    obs <- t(psi)%*%Y
    ehat <- expectile(obs, alpha)
    lab <- ifelse(obs > ehat, 1, -1)
    return(lab)
}

#FUNCTION: principalDirection main function. 
#find the principal direction by iteratively updating psi and the weights
prdir <- function(Y, alpha, reset.opt = "random", ini.opt = "mean", lab.ini = NA, iter.tol = 10, reset.tol = 50,nc=NA){
    p = dim(Y)[1]
    n = dim(Y)[2]
    #initiate some label vectors
    if(ini.opt == "mean"){
        lab = reinit(n, opt = "mean", Y = Y,nc=nc)
    }else{
        if(is.na(lab.ini)){
            lab = reinit(n, opt = "random",nc=nc)
        }else{
            lab = lab.ini
        }
    }
    change = TRUE
    iter = 1
    reset = 0
	conv=1
    while(change){
        change = FALSE
        iter = iter + 1
        psi <- psihat(Y, lab, alpha)        
        newlab <- expdir(Y, psi, alpha)
        differ <- sum(lab!=newlab)
        change <- ifelse(differ > 0, TRUE, FALSE)
        lab <- newlab        
#        print(paste("iter = ", iter, " differ = ", differ))
#        print(head(psi))
        #can get stuck in local min. In this case, initialize 
        #to a random label
        if(iter > iter.tol){
            lab <- reinit(n, opt = reset.opt, lab = lab, alpha = alpha, Y = Y, psi = psi)
            iter = 1
            reset = reset+1
		change=TRUE
            if(reset > reset.tol){
                print(paste("WARNING: exceed reset level for alpha = ", alpha))
                conv=0 #algorithm did not converge
		    change=FALSE
		    #return(list(psi,conv))
		}            
        }
    }
    print(paste("iter = ", iter, "reset = ", reset))
    return(list(psi,conv,differ,lab))
}

#FUNCTION: options for reset and initialize the labels
#random
#flip far away
#jiggle tau
reinit <- function(n, opt = "random", lab = NA, alpha = NA, Y = NA, psi = NA,nc=NA){
    if(opt == "random"){
        lab=ifelse(runif(n,0,1)>0.5, 1, -1)
    }
    if(opt == "far"){
        n = length(lab)
        lab=lab*ifelse(runif(n,0,1)>0.8, 1, -1)
    }
    if(opt == "jiggle"){
        delta = (0.5 - abs(alpha))/10
        alpha.new = alpha + runif(1,-delta, delta)
        lab=expdir(Y, psi, alpha.new)
    }
    if(opt == "mean"){
        C = cov(t(Y))
        pc1 = eigen(C)$vectors[,1]
        lab=expdir(Y, pc1, 0)
    }
	#check if at least one element of lab has the opposite sign, if not change sign of one elem
	if (all(lab==1) || all(lab==-1)){
		ind <-sample.int(n, size = nc, replace = FALSE)
		lab[ind]<-lab[ind]*(-1)
	}
	return(lab)
}


#FUNCTION: compute the first k principal components
pec.k <- function(Y, alpha, nk=2, reset.opt = "random", ini.opt = "mean", lab.ini = NA, iter.tol = 10, reset.tol = 50){
    p = dim(Y)[1]
    n = dim(Y)[2]
    psi.mat = matrix(NA, nrow = p, ncol = nk)
    mu.vec = rep(NA, nk)
    resid.obs <- Y
	conv=1
	cconv=NULL
    while(nk > 0){
	  print(paste("computing for k = ", nk))
        outp <- prdir(resid.obs, alpha, reset.opt = reset.opt, ini.opt = ini.opt, lab.ini = lab.ini, iter.tol = iter.tol, reset.tol = reset.tol,nc=nk)
        psi<-outp[[1]]; conv=outp[[2]]; lab=outp[[4]]; cconv=c(cconv,conv)
        psi.mat[,nk] <- psi
        obs = t(psi)%*%resid.obs
        mu.vec[nk] <- expectile(obs, alpha)
        resid.obs = resid.obs - as.matrix(psi)%*%obs	
	nk = nk-1
    }
    psi.mat = psi.mat[,rev(1:dim(psi.mat)[2])]
	return(list(rev(mu.vec), psi.mat, min(cconv),lab))
}

#FUNCTION: compute one component for a sequence of alpha
pec.alpha <- function(Y, alphaseq = seq(0.05, 0.45, by = 0.05), reset.opt = "random", ini.opt = "mean", iter.tol = 10, reset.tol = 50,nc=NA){
    p <- dim(Y)[1]
    n <- dim(Y)[2]
    psi.mat <- NA
    mu.vec <- NA
    lab.ini = NA
    for(alpha in alphaseq){
        print(paste("computing for alpha = ", alpha))
        psi = suppressWarnings(prdir(Y, alpha, reset.opt = reset.opt, ini.opt = ini.opt, lab.ini = lab.ini, iter.tol = iter.tol, reset.tol = reset.tol,nc=nc))
        psi.mat <- cbind(psi.mat, psi) 
        mu.vec <- append(mu.vec, expectile(t(psi)%*%Y, alpha))
        lab.ini <- expdir(Y, psi, alpha)
    }
    psi.mat <- psi.mat[,c(-1)]
    mu.vec <- mu.vec[-1]
    dimnames(psi.mat)[[2]] <- alphaseq
    return(list(mu.vec,psi.mat))
}

#FUNCTION: compute the residuals after projecting onto the subspace spanned by psi.mat
resid.proj <- function(Y, psi.mat,alpha){
        return(Y - psi.mat%*%t(psi.mat)%*%Y-apply(Y,1,mean))   
}

resid.pec <- function(Y,output.pec, tau, rss=FALSE, fitted=FALSE, mean=FALSE){
			psi.mat=as.matrix(output.pec[[2]])	
			mu.vec=output.pec[[1]]	
			d=dim(psi.mat)
			m=matrix(NA,d[1],d[2])
			labs=matrix(NA,ncol(Y),d[2])
			remain<-Y
			remain=remain-psi.mat%*%t(psi.mat)%*%remain
			rest=apply(remain,1,expectile,alpha=tau-0.5)
			remain=remain-rest
			if(rss){
				wij <-  sign(remain)*(tau - 0.5) + 0.5
				return(sum(remain^2*wij))
			}
	        	if(fitted){
				return(Y-remain)
			}
			if(mean){
				return(apply(Y,1,mean)+rest)
			}

			else{
				return(remain)
			}
   
}


#FUNCTION: compute pec for a sequence of alpha and k components?
pec.main <- function(Y, alphaseq = seq(0.05, 0.45, by = 0.05), k = 1, reset.opt = "random", ini.opt = "mean", iter.tol = 10, reset.tol = 50){
    p <- dim(Y)[1]
    n <- dim(Y)[2]
    a <- length(alphaseq)
    psi.array <- array(NA, c(p, k, a))
    mu.mat <- array(NA, c(k, a))
    lab.ini = NA
    for(i in c(1:a)){
        alpha = alphaseq[i]
        print(paste("computing for alpha = ", alpha))
        out = suppressWarnings(pec.k(Y, alpha = alpha, nk = k, reset.opt = reset.opt, ini.opt = ini.opt, lab.ini = lab.ini, iter.tol = iter.tol, reset.tol = reset.tol))
        psi.array[,,i] <- out[[2]]
        mu.mat[,i] <- out[[1]]
        lab.ini <- expdir(Y, psi.array[,1,i], alpha)
    }
    return(list(mu.mat, psi.array))
}

#FUNCTION: plot a sequence of vectors on the same plot
#the vectors are columns of the supplied matrix

plot.vec <- function(mat, main.str = "",type="l",pch=1,add=FALSE,ylim=NA,lwd=2,color=rainbow,colp="gray80",xlab="",ylab="")
{
	a=1
	if(any(is.na(ylim))){ylim=c(min(mat,na.rm=T),max(mat,na.rm=T))}
	if (add==FALSE){ 
    		a=2	
		plot(mat[,1], type = type, lwd=lwd, ylim = ylim, main = main.str, pch=pch,col=colp,axes=FALSE,xlab=xlab,ylab=ylab)
    	}
	for(i in c(a:dim(mat)[2])){
		if (type=="p"){
		points(mat[,i], col = colp,pch=pch)
		} else { lines(mat[,i], col = color(i),lwd=lwd)}
    } 
}

#define gray color
expgr<-function(n){
col="gray80"
return(col)
}


