
setwd("C:/Users/Vitória/Dropbox/FunctionalSeriesPackage/code/newcode/codefinal2")
source('newcode_para_rodar.r')


N_repl = 100 # COLOCA = 100 OU = 1000
sigma_grid.length = 11 # COLOCA = 11

heatmap_array = array(0, c(sigma_grid.length,sigma_grid.length,N_repl))


for (repl in 1:N_repl){
	for (p in c(5)){
		for (q in c(.5)){
			for (n in c(100)){
				#heatmap_matrix=matrix(NA,nrow = length(seq(from=1,to=2,length=11)), ncol = length(seq(from=1,to=2,length=11)))
				i=0
				for(sigma_loop1 in seq(from=1,to=2,length=sigma_grid.length)){
					j=0
					i =i+1
					for(sigma_loop2 in seq(from=.5,to=1,length=sigma_grid.length)){
						j =j+1
						pvalue1 = pvalue2 = .5
						c=0
						while(pvalue1 + pvalue2 > .05){
						  c = c+1
							resultados=mcsimu(n=n, sigma1=sigma_loop1,sigma2=sigma_loop2,q =q,p = p)
              pvalue1 = resultados$pvalue1
              pvalue2 = resultados$pvalue2
							if (c>1000) break()
						}
						#heatmap_matrix[i,j]=resultado_final
						heatmap_array[i,j,repl] = resultados$correct_classification_rate
						#cat('\nrepl = ', repl, '\np = ', p, '\nq = ', q, '\nn = ', n, '\nsigma1 = ', sigma_loop1, '\nsigma2 = ', sigma_loop2,
						#		file = 'foo.txt', append = FALSE)
					}
				}
			}
		}
	}
}



heatmap_matrix = matrix(0, sigma_grid.length, sigma_grid.length)
for( i in 1:sigma_grid.length){
	for( j in 1:sigma_grid.length){
		heatmap_matrix[i,j] = mean(heatmap_array[i,j,])
	}
}


library(plotly)
plot_ly(z = heatmap_matrix, type = "heatmap")
