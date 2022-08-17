super_fun = function(u, Y, dimension = 11, p = 1){
	# Setting up the functions
	# Defining the integral on L^2([a,b])
	integral = function(f, u, a, b) integrate(splinefun(u, f), a, b)$value
	
	# Defining the inner product on L^2([a,b])
	inner.product = function(f, g, dmu) drop(crossprod(f, g))*dmu
	
	# Defining the L2-norm
	L2norm = function(f, dmu) sqrt(inner.product(f, f, dmu))
	
	# Defining a function sampleCols, which samples the columns of any given matrix
	sampleCols = function(A){
		idx = sample(1:ncol(A))
		M = matrix(0, nrow = nrow(A), ncol = ncol(A))
		for (i in 1:length(idx)) {
			M[,i] = A[,idx[i]]
		}
		M = M
	}
	
	du = u[2]-u[1]
	m = nrow(Y)
  n = ncol(Y)
# d0 is the hypothesized dimension
## RUN THE BOOTSTRAP IN ORDER TO SELECT d0 PROPERLY
d0 = dimension

# Defining Ybar
  Ybar = rowMeans(Y)

# Defining the deviation function Ydev = Y - Ybar, which is used as an input
# by the function inner.product in constructing the matrix Kstar
  Ydev = Y - Ybar

  core = inner.product(Ydev,Ydev, du)
  Kstar.core0 = core[1:(n-p),1:(n-p)]
  Kstar.core = array(0,c(n-p,n-p,p))
  for (k in 1:p) Kstar.core[,,k] = core[(k+1):(n-(p-k)),(k+1):(n-(p-k))]

# Summing the matrices in 'Kstar.core'
  Kstar.sum = matrix(0,nrow=n-p,ncol=n-p)
  for (k in 1:p) Kstar.sum = Kstar.sum + Kstar.core[,,k]

# Defining Kstar
  Kstar = (n-p)^(-2) * Kstar.sum %*% Kstar.core0

# 4. Eigen-analisys

# Getting the eigenvalues and eigenvectors from 'Kstar'

# Storing the eigenvalues; length(thetahat)=n-p
  thetahat = eigen(Kstar)$values

# Storing the eigenvectors; each column of gammahat corresponds to
# one eigenvector; dim(gammahat)=[(n-p) x (n-p)]
  gammahat = eigen(Kstar)$vectors

# Defining a tolerance level for 'testing' if the imaginary part of the
# eigenvalues and eigenvectors is zero
tol = 10^(-4)

# Checking if there are any complex eigenvalues (among the first eleven)
for (j in 1:11){
  if (abs(Im(thetahat[j]))>tol) print("Complex eigenvalue found.")
}
  
  thetahat = Re(thetahat)
  thetahat = sort(thetahat, index.return=TRUE, decreasing=TRUE)
  thetahat.index = thetahat$ix
  thetahat = thetahat$x

# Ordering the eigenvectors accordingly
  gammahat.temp = matrix(0,nrow=nrow(gammahat),ncol=ncol(gammahat))
  for (j in 1:(length(thetahat))){
    gammahat.temp[,j] = gammahat[,thetahat.index[j]]
  }
  gammahat = gammahat.temp

# Storing the original eigenvalues and eigenvectors
  thetahat.old = thetahat
  gammahat.old = gammahat

# 5. Defining the estimator Yhat

# Storing only the d0 largest eigenvalues and the associated eigenvectors
  thetahat = Re(thetahat.old[1:d0]) # length(thetahat) = d0
  gammahat = Re(gammahat.old[,1:d0]) # dim(gammahat) = [(n-p) x d0]

# Defining the eigenfunctions psihat.root. These are the functions
# given in equation (2.13)
#   psihat.root[u,j] = gammahat[1,j]*Ydev[u,1] + ... + gammahat[n-p,j]*Ydev[u,n-p]
# Note that dim(psihat.root) = [m x d0]. These functions are not
# necessarily orthonormal.
  psihat.root = Ydev[,1:(n-p)]%*%gammahat

# Normalizing (tested for orthogonality; already orthogonal to each other)
  psihat = matrix(0,nrow=m,ncol=d0)
  #psihat = matrix(0,nrow=length(alpha_grid),ncol=d0)
  for (i in 1:d0) psihat[,i] = psihat.root[,i]/L2norm(psihat.root[,i], du)

# Defining etahat (d0 x n). We have that
#   etahat[j,t]=<Ydev[,t],psihat[,j]>.
# This is a d0-vector process. Time varies columnwise.
  etahat = inner.product(psihat,Ydev, du)

# Defining Yhat (m x n). The (u,t)-th element of Yhat is given by
#   Yhat[u,t] = Ybar[u] + psihat[u,1]*etahat[1,t] + ... + psihat[u,d0]*etahat[d0,t]
  Yhat = Ybar + psihat%*%etahat

# # Note that no restrictions were made upon Yhat, so it may happen that
# sum(Yhat[,t])*du!=1 and/or Yhat[u,t]<0.
# We define Yhat.fix to meet these restrictions.
Yhat.fix = Yhat
Yhat.fix[Yhat.fix<0] = 0
for (t in 1:n) Yhat.fix[,t] = Yhat.fix[,t]/(sum(Yhat.fix[,t])*du); rm(t)

# Defining epsilonhat (m x n)
  epsilonhat = Y - Yhat
  return(
    list(
      Y = Y,
      Ybar = Ybar,
      thetahat = thetahat,
      gammahat = gammahat,
      psihat = psihat,
      etahat = etahat,
      Yhat = Yhat,
      Yhat.fix = Yhat.fix,
      epsilonhat = epsilonhat,
      u = u
    )
  )
}

#fooY = super_fun(Y, 11)
