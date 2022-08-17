# COMENTÁRIOS MARCADOS COM DUPLA INTERROGAÇÃO '??'
#Performs Bathias model and return the dim_max graphs for the auto-functions and 1 graph for the
#mean of the estimate of y

Ftsa2=function(data=t(Y),  endpoints=c(18.0, 29.2),  n_grid=ncol(data),  p=6,  dim_max=2,  density_fix=TRUE
  ,Shiny=FALSE){
  
      if(is.matrix(data)==FALSE)stop("Your object data must be of matrix type.") 
      if(p<1)stop("p must be greater than 0.")
      if(n_grid<1)stop("n_grid must be greater than 0.") 
      if(dim_max<1)stop("dim_max must be greater than 0.")
  
  a=endpoints[1]
  b=endpoints[2]
  
      if(b<=a)stop("b must be greater than a.") 

      data=t(data)
      Kbar=rowMeans(data)
      Kdev = data - Kbar
      du=(b+abs(a))/n_grid

  inner_product=drop(crossprod(Kdev))*du

  m=ncol(inner_product)
  Kstar.core0 = inner_product[1:(m-p),1:(m-p)]
  Kstar.sum = matrix(0,nrow=m-p,ncol=m-p)
  for (k in 1:p) {
    Kstar.sum = Kstar.sum + inner_product[(k+1):(m+k-p),(k+1):(m+k-p)]
  }

  Kstar = (m-p)^(-2) * crossprod(Kstar.sum, Kstar.core0)

      thetahat = as.vector(eigen(Kstar)$values)
      gammahat = as.matrix(eigen(Kstar)$vectors)
      thetahat = Re(thetahat)
      thetahat = sort(thetahat, index.return=TRUE, decreasing=TRUE)
      thetahat.index = as.vector(thetahat$ix)
      thetahat = as.vector(thetahat$x)


  gammahat.temp = matrix(0,nrow=nrow(gammahat),ncol=ncol(gammahat))
  for (j in 1:(length(thetahat))){
    gammahat.temp[,j] = gammahat[,thetahat.index[j]]
  }

      thetahat.old = thetahat
      gammahat.old = gammahat
      thetahat = Re(thetahat.old[1:dim_max])
      gammahat = Re(gammahat.old[,1:dim_max])

  psihat.root = drop(crossprod(t(Kdev[,1:(m-p)]),gammahat))

  psihat = matrix(0,nrow=n_grid,ncol=dim_max)
  for (i in 1:dim_max) {
    psihat[,i] = psihat.root[,i]/sqrt(drop(crossprod(psihat.root[,i]))*du)
  }

      etahat = crossprod(Kdev,psihat)*du

      Yhat = Kbar + crossprod(t(psihat),t(etahat))

  Yhat.fix = Yhat
  
  if (density_fix==TRUE) Yhat.fix[Yhat.fix<0] = 0

  for (t in 1:m) {
    Yhat.fix[,t] = Yhat.fix[,t]/(sum(Yhat.fix[,t])*du)
  }

      epsilonhat = data - Yhat.fix

  for (i in 1:length(psihat[,1])){ # o que esse loop faz?
    j=length(psihat[,1])+1-i
    if(abs(psihat[i,1])>10*exp(-4)) end=i
    if(abs(psihat[j,1])>10*exp(-4)) begin=j
  }
      
      if(a<0 && b<0){
        end=length(psihat[,1])
        begin=1}
      if(a>0 && b>0){
        end=length(psihat[,1])
        begin=1}

      scale_matrix=matrix(0,ncol=ncol(psihat)+1,nrow = end-begin+1)
      scale_matrix[,1]=Kbar[begin:end]
      scale_matrix[,2:ncol(scale_matrix)]=psihat[begin:end,]
      scale_matrix_new=scale(scale_matrix,center = FALSE,scale = TRUE)

      min_vector=vector()
      max_vector=vector()

  for (k in 1:ncol(scale_matrix_new)){
    min_vector[k]=min(scale_matrix_new[,k])
    max_vector[k]=max(scale_matrix_new[,k])
  }

      Plot_lower_bound=min(min_vector)
      Plot_upper_bound=max(max_vector)


  if (ncol(psihat)>0){

    layout(matrix(c(1,2,3,4,5,6), ncol = 3, byrow = TRUE))

    x = seq(a,b-du,du)[begin:end]
    Ybar=Kbar[begin:end]

  }

  if (Shiny==TRUE){
    
  }    
      
      
  lista=list(data,thetahat.old,gammahat.old,epsilonhat,Kdev,Kbar,etahat,psihat)
  return(lista)
}
