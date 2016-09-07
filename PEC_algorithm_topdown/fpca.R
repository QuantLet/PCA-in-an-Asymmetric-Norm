# ------------------------------------------------------------------------------
# Book:       	N.M. Tran, M. Osipenko, and W.K. Härdle 
#               "Principal Components in an Asymmetric Norm"
# ------------------------------------------------------------------------------
# Quantlet:   	fpca.R
# ------------------------------------------------------------------------------
# Description:	Functional library for the TopDown and BottomUp algorithms
# ------------------------------------------------------------------------------
# See also:     China_example.R, PrincipalExpectile.R
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



#--Functional library for the TopDown and BottomUp algorithms

#--- laws internal one-step function
laws.nxt <- function
###Internal function. Given a basis, compute the next component in LAWS
###Return the UN-NORMALIZED coefficients v and weights wij
(X, 
###data. p times n matrix, p is the data dimension, n the number of observations
B, 
###basis. p times k matrix.
tau, 
###a number in (0,1)
wij.ini, 
###initial weight matrix. p times n matrix.
mu = NA, 
###We center X using mu, if supplied.
max.iter = 30, 
###max number of iterations allowed to try before converging
messages = FALSE
###If true, we will output debugging messages
){
	p <- nrow(X)
	n <- ncol(X)
	k <- ncol(B)
	#center X using the given "mean" vector
	if(!is.na(mu[1])){
		X <- apply(X, 2, '-', mu)
	}
	diff <- 1
	iter <- 1
	wij <- wij.ini
	while(diff > 0 && iter < max.iter){
		#compute for each vi: vi = (B^tWiB)^{-1}B^tWiYi
		v <- matrix(NA, nrow = k, ncol = n)
		width <- length(wij[,1])
		for(j in c(1:n)){
			W <- diag(wij[,j], nrow = width, ncol = width)
			v[,j] <- solve(t(B)%*%W%*%B)%*%t(B)%*%W%*%X[,j]
		}
		#update the weights
		resid <- X - B%*%v
		wij.new <- sign(resid)*(tau - 0.5) + 0.5 
		diff <- sum(wij != wij.new)
		#copy over the updated wij
		wij <- wij.new
		#update iteration number
		iter = iter +1
	}
	
	if(messages){
		if(diff == 0){
			print(paste("converged in ", iter, "steps"))
		}
		else{
			print(paste("did not converge in ", max.iter, "steps"))
		}
	}
	return(list(v, wij))
}


laws.const <- function
###Public function. Estimates the constant term, that is, the tau-expectile.
(X, 
###data. p times n matrix, p is the data dimension, n the number of observations
tau, 
###A number in (0,1)
max.iter = 30, 
###Max number of iterations allowed
tol = 1e-10
###Error tolerance level
){
	p <- nrow(X)
	n <- ncol(X)
	basis <- as.matrix(rep(1, n), ncol = 1)
	#set initial weights to be constant
	wij <- matrix(0.5, nrow = p, ncol = n)	
	#find the constant term
	out <- laws.nxt(t(X), basis, tau, t(wij), mu = NA, max.iter, messages = FALSE)
	mu <- out[[1]]
	#return mu, the constant term approximation
	return(mu)
}

#------------ main laws function
laws.main <- function
###Public function. Uses LAWS to compute a low-rank approximation that minimizes the ell_{2,tau}-error.
(Y, 
###data. p times n matrix, p is the data dimension, n the number of observations
nc, 
###Number of components to compute. Also is the rank of the low-rank approximation.
tau, 
###A number in (0,1) 
max.iter = 30, 
###Max number of iterations allowed
max.reset = 50, 
###Max number of resets of initial values allowed
tol = 1e-10, 
###Error tolerance level
mu = NA, 
###We center Y using mu. If not pre-supplied, the tau-expectile will be used.
mu.estimate = TRUE,
###If TRUE, the constant term will be re-estimated
preB = NA,
###basis vectors pre-supplied, required to have norm 1. If is not NA, this will forces the low-rank approximation to contain this basis. In this case, nc will be the number of extra components. The total rank of the approximation will be rank(preB) + nc.
basis.start = NA
###basis vectors to initiate the algorithm, if NA chosen to be random columns of X
){	
	ind=1#indicator of convergence
	p <- nrow(Y)
	n <- ncol(Y)
	if(!any(is.na(preB))){#changed to omit the warning
		prek <- dim(preB)[2]
		bsupplied <- TRUE
	} else{
		prek <- 1
		bsupplied <- FALSE
	}
	k <- nc+prek+1
	if(is.na(mu[1])){
		mu <- laws.const(Y, tau)
	}
	#basis to initiate
	if(!any(is.na(basis.start))){#changed to omit the warning
		basis <- basis.start
	} else{
		#initial extra basis is chosen to be RANDOM nc columns of Y, normalized
		col <- sample.int(n, size = nc, replace = FALSE)
		basis <- scale(Y[,col], center=FALSE)/sqrt(p-1) #basis is p times nc
	}
	#set initial weights to be constant
	wij<- matrix(0.5, nrow = p, ncol = n)	
	err<-1
	iter<-1
	reset<-0
	while(err > tol){
		#first iteration: find the vi and the weight wij
		if(!bsupplied){
			out <- laws.nxt(Y, B = basis, tau = tau, wij = wij, mu = mu, max.iter, messages = FALSE)
			v <- t(out[[1]]) #v is n times k-1
			wij <- out[[2]]
			#second iteration: fix [1, v] to be the basis (don't need to normalize v)
			if(mu.estimate){ #if need to estimate mu
				v <- cbind(rep(1, n), v)
			}
			out <- laws.nxt(t(Y), B = v, tau = tau, wij = t(wij))
		} else{
			out <- laws.nxt(Y, cbind(preB,basis), tau, wij, mu = mu, max.iter, messages = FALSE)
			wij <- out[[2]]
			#subtract off the coefficients of preB
			vprek <- t(out[[1]])[,1:prek]
			Yprime <- Y - preB%*%t(vprek)
			v <- t(out[[1]])[,(prek+1):(k-1)]
			if(mu.estimate){
				v <- cbind(rep(1, n), v)
			}
			out <- laws.nxt(t(Yprime), v, tau, t(wij))
		}
		if(mu.estimate){ #if need to estimate mu
			mu <- t(out[[1]])[,1]	
			basis.new <- t(out[[1]])[,2:(nc+1)]
		} else{
			basis.new <- t(out[[1]])[,1:nc]
		}
		#update wij
		wij <- t(out[[2]])		
		#normalize the basis vectors
		basis.new <- scale(basis.new, center=FALSE)/sqrt(p-1)
		#check for convergence of the basis vector, 
		#taking into account the potential flip in sign
		err <- min(sum((basis.new - basis)^2), sum((basis.new + basis)^2))
		basis <- basis.new	#update the basis
		iter <- iter + 1
		if(iter > max.iter){
			col <- sample.int(n, size = nc, replace = FALSE)
			basis <- scale(Y[,col], center=FALSE)/sqrt(p-1) #basis is p times nc
			#reset initial weights to be constant
			wij<- ((tau - 0.5) + 0.5)*matrix(ifelse(runif(n*p,0,1)>0.5,1,-1),nrow = p, ncol = n)#matrix(0.5, nrow = p, ncol = n)
			iter=1
			reset=reset+1
			if (reset > max.reset){
				print(paste("WARNING: exceed reset level for tau = ", tau))
				break
			}
		}
	}
	if(err > tol){
		print(paste("Warning: did not convergeã€€in ", max.iter, " iterations"))
		ind=0 #no convergence
		#stop("no convergence")
		
	}
	#return the main stuff
	if(any(is.na(preB))){ #changed to omit the warning
		return(list(v, mu, basis, wij,ind))	
	} else{
		vfinal <- cbind(vprek, v[,2:(nc+1)])
		vfinal <- cbind(rep(1,n), vfinal)
		return(list(vfinal, mu, cbind(preB, basis), wij, ind))
	}
}



#FUNCTION: BottomUpFast
pp <- function
###public function BottomUpFast. Projection pursuit up to the first k components.
###Returns a list in the following order: coefficients, constant term mu, basis vector, weight matrix wij
(X,
###data. p times n matrix, p is the data dimension, n the number of observations
k, 
###Number of components to compute.
tau 
###A number in (0,1)
){
	output <- laws.main(X, nc = 1, tau = tau, max.iter = max.iter, tol = tol)
	coeff <- output[[1]]
	mu <- output[[2]]
	basis <- output[[3]]	
	for(i in 2:k){
		X <- resid(X, output) #compute residual
		#compute next decomposition
		output <- laws.main(X, nc = 1, tau = tau, max.iter = max.iter, tol = tol)
		coeff <- cbind(coeff, output[[1]][,2])
		mu <- mu + output[[2]] #modified to estimate the constant as in pp()
		basis <- cbind(basis, output[[3]])	
	}
	X <- resid(X, output)
	wij <- sign(X)*(tau - 0.5) + 0.5
	#normalize the basis vectors
	if(k > 1){	
        coeff <- coeff[,2:(k+1)]%*%t(qr.R(qr(basis)))
        basis <- qr.Q(qr(basis))
    }
	return(list(coeff, mu, basis, wij))
}

#FUNCTION: BottomUp
bup <- function
###public function BottomUp. Run projection pursuit, but recompute the coefficients of the previously found basis when adding on the next basis vector.
###Returns a list in the following order: coefficients, constant term mu, basis vector, weight matrix wij
(Y,
###data. p times n matrix, p is the data dimension, n the number of observations
k, 
###Number of components to compute.
tau,
###A number in (0,1)
max.reset=50,
basis.start = NA
){
	output <- laws.main(Y, nc = 1, tau = tau,max.reset=max.reset, basis.start=basis.start)
	basis <- output[[3]]
	if(k > 1){	 #added if, otherwise computes 3 PC for k=1
		i=2
		while(i <= k && output[[5]]==1){
			output <- laws.main(Y, nc = 1, tau = tau, preB = basis, max.reset=max.reset, basis.start=basis.start)
			basis <- output[[3]]
			i=i+1
		}
		
	  #normalize the basis vectors
	  #output[[1]] <- output[[1]][,2:i]%*%t(qr.R(qr(output[[3]])))
        output[[1]][,2:i] <- output[[1]][,2:i]%*%t(qr.R(qr(output[[3]]))) #CHANGED HERE to keep column of ones for mean -> fun resid(); i=k+1
        output[[3]] <- qr.Q(qr(output[[3]]))

    }
	return(output)
}

#--- FUNCTION: residual statistics
resid <- function
###Returns statistics on the residuals from a laws.main output. Mainly used for debugging.
###The default option returns the residual vector.
(X, 
###data. p times n matrix, p is the data dimension, n the number of observations
output, 
###A list obtained as the output of bup or pp
tau,
###the expectile level
rss = FALSE, 
###If true, returns the residual sum of squares.
fitted = FALSE
###If true, returns the fitted values.
){
	coeff <- output[[1]]
	basis <- cbind(output[[2]], output[[3]])	
	remain <- X - basis%*%t(coeff)
	if(rss){
		wij <-  sign(remain)*(tau - 0.5) + 0.5
		return(sum(remain^2*wij))
	}
	if(fitted){
		return(X-remain)
	}
	else{
		return(remain)
	}
}

#---- FUNCTION: compare basis
#This is to compare the TopDown output with BottomUp and PrincipalExpectile
compareBasis <- function
### Project basis.old onto the span of a fixed basis basis.fix.
(basis.fix, 
###Fixed basis
basis.old
###Basis to be re-expressed in the span of basis.fix
){
    qrobj <- qr.solve(basis.fix, basis.old)
    basis.new <- basis.old%*%solve(qrobj)
    return(basis.new)
}

compareBasisVec <- function
### Project vectorized version of basis.old onto the span of a fixed basis basis.fix.
(basis.fix, 
###Fixed basis
basis.old.vec,
###Basis to be re-expressed in the span of basis.fix
nc=2 
###Number of dimensions
){	
	p=length(basis.old.vec)/nc
	basis.old <- matrix(basis.old.vec,p,nc)
	qrobj <- qr.solve(basis.fix, basis.old)
	basis.new <- basis.old%*%solve(qrobj)
	return(matrix(basis.new,nc*p,1))
}

### function to compare bases returned in a list
listBasis<-function(basis0, listb, num){
###basis0 is the fixed basis, listb is a list with outputs from TD, BUP or PEC, num is the location of basis in the output (num=3 for BUP,TD, and num=2 for PEC)
	return(as.matrix(compareBasis(basis0,listb[[num]])))
}



