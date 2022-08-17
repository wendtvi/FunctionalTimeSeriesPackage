
#setwd("C:/Users/Pichau/Dropbox/FunctionalSeriesPackage/code/newcode")

source('newcode_para_rodar.r')


N_repl = 10 # COLOCA = 100 OU = 1000
sigma_grid.length = 3 # COLOCA = 11

heatmap_array = array(0, c(sigma_grid.length,sigma_grid.length,N_repl))


for (repl in 1:N_repl){
  for (p in c(5)){
    for (q in c(.5)){
      for (n in c(100)){
        heatmap_matrix=matrix(NA,nrow = length(seq(from=1,to=2,length=11)), ncol = length(seq(from=1,to=2,length=11)))
        i=0
        for(sigma_loop1 in seq(from=1,to=2,length=sigma_grid.length)){
          j=0
          i =i+1
          for(sigma_loop2 in seq(from=.5,to=1,length=sigma_grid.length)){
            j =j+1
            resultados=mcsimu(n=n, sigma1=sigma_loop1,sigma2=sigma_loop2,q =q,p = p)
            c=0
            resultado_teste=mean(resultados$thetahat.fts)
            while(mean(resultados$thetahat.fts)<=.65){
              c=c+1
              resultados=mcsimu(n=n, sigma1=sigma_loop1,sigma2=sigma_loop2,q =q,p = p)
              if(mean(resultados$thetahat.fts)>resultado_teste) {
                resultado_final=mean(resultados$thetahat.fts)
              }
              if(mean(resultados$thetahat.fts)<=resultado_teste) resultado_final=resultado_teste
              if (c==30) break()
              resultado_teste=resultado_final
            }
            #heatmap_matrix[i,j]=resultado_final
            heatmap_array[i,j,repl]
          }
        }
        str_name=paste("p", p, "q", q, "n", n)
        write.csv2(heatmap_matrix,file=str_name)
      }
    }
  }
}


