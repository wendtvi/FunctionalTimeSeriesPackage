library(depmixS4)
# setwd("C:/Users/Pichau/Dropbox/FunctionalSeriesPackage/code/code-vitoria/functions")
source('4 Ftsa2 - comments.R')
cpvalue1=0
cpvalue2=0

mcsimu = function(n=200, # sample size
									sigma1=.8, # standard deviation of first error term # ?? VARIAR DE 1 a 2 NUMA GRADE COM 11 PONTOS seq(from=1,to=2,length=11)
									sigma2=1.2, # standard deviation of second error term # ?? VARIAR DE 0.5 a 1 NUMA GRADE COM 11 PONTOS seq(from=.5,to=1,length=11)
									q = .5, # transition probability
									p = 20 # lag-covariance parameter (see Bathia et al. (2010))
					        ){
  
	# PARTE 1: EXPLORANDO O MODELO (1 UNICA REALIZACAO, GRAFICOS, ETC)
	burn_in = 1000
	STATES = c(0,1) # ESPACO DE ESTADOS DA CADEIA DE MARKOV LATENTE
	
	Q = rbind(c(1 - q, q),
						c(q, 1 - q))
	
	invariant_dist = c(1/2, 1/2)
	
	# IMPORTANTE: RENOMEEI A VARIAVEL LATENTE PARA THETA, ASSIM FICA MENOS CONFUSO (SE NAO, O LEITOR PODERIA PENSAR QUE ERA SEMPRE A MEDIA)
	theta0 = sample(STATES, size = 1, prob = invariant_dist)
	
	theta = numeric()
	
	# IMPORTANTE: ALTEREI AS DIMENSOES DE Z: AGORA O TEMPO t PERCORRE AS COLUNAS. ISSO FICA COERENTE COM A SAIDA DO ftsa2
	Z = matrix(0, nrow = 2, ncol = n+burn_in)
	
	mu = function(theta){
		c(4*theta-2, 2*theta-1)
	}
	
	current_state = theta0
	
	for (t in 1:(n+burn_in)){
		theta[t] = sample(STATES, size = 1, prob = Q[current_state+1,])
		Z[,t] = mu(theta[t]) + rnorm(2,0,sd = c(sigma1, sigma2))
		current_state = theta[t]
	}
	
	z1 = Z[1, (burn_in+1):(burn_in+n)]
	z1.lag = Z[1, burn_in:(burn_in+n-1)]
	z2 = Z[2, (burn_in+1):(burn_in+n)]
	z2.lag = Z[2, burn_in:(burn_in+n-1)]
	
	
	theta = theta[(burn_in+1):(n+burn_in)]

# ?? Essa parte podemos pular: vamos fazer a classificacao usano apenas os loadings estimados
	# df = as.data.frame(list(z1,z1.lag,z2,z2.lag))
	# fmsp = depmix(list(z1~z1.lag+z2.lag, z2~z1.lag+z2.lag),nstates=2,data=df,family=list(gaussian(),gaussian()))
	# 
	# foo = try(fit(fmsp))
	# 
	# thetahat=vector()
	# test1=posterior(foo)[,1] - 1
	# test2 = abs(test1-1)
	# #test2=replace(posterior(foo)[,1],posterior(foo)[,1]==1,1)
	# #test2=replace(posterior(foo)[,1],posterior(foo)[,1]==2,0)
	# 
	# eqm1=mean((theta-test1)^2)
	# eqm2=mean((theta-test2)^2)
	# 
	# if(eqm1>=eqm2) thetahat=test2
	# if(eqm1<eqm2) thetahat=test1
	
	
	# GERANDO A SERIE FUNCIONAL
	# 
	# AUTOFUNCOES
	# 
	
	# PRIMEIRA AUTOFUNCAO
	phi1 = function(u) sqrt(2)*sin(2*pi*u)
	
	# SEGUNDA AUTOFUNCAO
	phi2 = function(u) sqrt(2)*cos(2*pi*u)
	
	X = function(u,t)  z1[t]*phi1(u) + z2[t]*phi2(u)
	
	
	# GRADE DE PONTOS NO DOMINIO DAS FUNCOES
	du = .01
	u = seq(from=0,to=1,by=du)
	
	# AVALIANDO AS FUNCOES NA GRADE DE PONTOS
	phi1_eval = phi1(u)
	phi2_eval = phi2(u)
	X_eval = sapply(1:n, function(t){
		X(u, t)
	})
	
	
	# ADICIONANDO RUIDO (QUE NESSE EXEMPLO E' UMA BROWNIAN BRIDGE)
	BrownianBridge = sapply(1:n, function(t){
		foo = rnorm(length(u), sd = du)
		foo = cumsum(foo)
		return(foo - u*foo[length(u)])
	})
	
	scaleBB = 1 # POR ENQUANTO NAO VAMOS MEXER NESSE PARAMETRO
	BrownianBridge = BrownianBridge*scaleBB
	
	# SERIE TEMPORAL FUNCIONAL EFETIVAMENTE OBSERVADA:
	Y_eval = X_eval + BrownianBridge
	
	# APLICANDO A METODOLOGIA DO BATHIA ET AL (2010)
	# POR ENQUANTO USEI O MEU CODIO ANTIGO PORQUE ESTOU MAIS HABITUADO AS PARAMETRIZACOES E AS SAIDAS DA FUNCAO. FICA A VONTADE PRA MUDAR
	foon = Ftsa2(t(Y_eval), endpoints=c(0, 1),  n_grid=ncol(t(Y_eval)),  p=p,  dim_max=2,  density_fix=FALSE
	             ,Shiny=FALSE)
	
	etahat=foon[[7]]
	etahat=as.matrix(etahat)
	# E' PRECISO FAZER ALGUNS AJUSTES PARA VISUALIZAR AS AUTOFUNCOES E AS SERIES TEMPORAIS DOS COEFICIENTES ESTIMADAS (ISSO VARIA DE SEED PRA SEED; NA SIMULACAO DE MONTE CARLO ABAIXO ISSO E' IRRELEVANTE POIS SO' VAMOS OLHAR PARA AS ESTIMATIVAS DA CADEIA DE MARKOV LATENTE)
	
	z1hat = etahat[,1]
	z1hat.lag = c(0,etahat[1:(n-1),1]) # acrscento um zero na primeira posicao da serie defasada
	z2hat = etahat[,2]
	z2hat.lag = c(0,etahat[1:(n-1),2])
	
	# APLICO DEPMIX NA SERIE TEMPORAL DOS COEFICIENTES
	df = as.data.frame(list(z1 = z1hat, z1.lag = z1hat.lag, z2 = z2hat, z2.lag = z2hat.lag))
	fmsp = depmix(list(z1~z1.lag+z2.lag, z2~z1.lag+z2.lag),nstates=2,data=df,family=list(gaussian(),gaussian()))
	
	
	m=-1
	while(m==-1){
	foo = try(fit(fmsp, emcontrol = em.control(maxit=1000, tol = 1e-17, crit = 'relative', random.start=TRUE, classification = 'soft')))
	if(typeof(foo)=="character"){ 
	  m=-1
	}
	if(typeof(foo)!="character"){ 
	thetahat.fts=vector()
  	test1=posterior(foo)[,1] - 1
  	# ?? aqui possivelmente esta' causando um erro. vamos tentar algo mais simples:
  	test2 = abs(test1-1)
  	#test2=replace(posterior(foo)[,1],posterior(foo)[,1]==1,1)
  	#test2=replace(posterior(foo)[,1],posterior(foo)[,1]==2,0)
  	
  	eqm1=mean((theta-test1)^2)
  	eqm2=mean((theta-test2)^2)
  	
  	if(eqm1>=eqm2) thetahat.fts=test2
  	if(eqm1<eqm2) thetahat.fts=test1
  	
  	m=Inf
  	

	}
	}
	
	
	pvalue1 = summary(lm(z1hat ~ thetahat.fts))$coef[2,4]
	pvalue2 = summary(lm(z2hat ~ thetahat.fts))$coef[2,4]
	
	if (typeof(try(summary(lm(z1hat ~ thetahat.fts))$coef[2,4]))!="character") {
	  pvalue1=0
	  cpvalue1=cpvalue1+1
	}
	if (typeof(try(summary(lm(z2hat ~ thetahat.fts))$coef[2,4]))!="character") {
	  pvalue2=0
	  cpvalue2=cpvalue2+1
	}
	
	return(
		list(
			theta = theta,
#			thetahat = thetahat,
			thetahat.fts = thetahat.fts,
			z1 = z1,
			z2 = z2,
			z1hat = z1hat,
			z2hat = z2hat,
			phi1 = phi1_eval,
			phi2 = phi2_eval,
			X = X_eval,
			Y = Y_eval,
			Bridge = BrownianBridge,
			u = u,
			du = du,
      correct_classification_rate = mean(abs(thetahat.fts - theta)<1/2),
      pvalue1 = summary(lm(z1hat ~ thetahat.fts))$coef[2,4],
      pvalue2 = summary(lm(z2hat ~ thetahat.fts))$coef[2,4]
		)
	)
}

	
	# As linhas de codigo abaixo utilizam outro pacote para fazer a classificacao por Regime Switching VAR
	# Ao que tudo indica, algumas replicacoes realmente apresentam maior dificuldade de classificacao
	# 
	# data=array(t(foon$etahat),dim=c(n,1,2))
	# theta.init = init.theta.MSAR(data,M=2,order=0)
	# mod.hh = fit.MSAR(data,theta.init,verbose=TRUE,MaxIter=10)
	# regimes=regimes.plot.MSAR(mod.hh,data,ylab="temperatures")
	# regimes[regimes<=1]=0
	# regimes[regimes>1]=1
	# 
	# plot(z1,z2, pch = 20, type = 'n')
	# points(z1[theta<1], z2[theta<1], pch = 19, col = rgb(.7,0,.2,.5),cex=1.4)
	# points(z1[theta>0], z2[theta>0], pch = 19, col = rgb(.2,0,.7,.5),cex=1.4)
	# 
	# points(z1[regi