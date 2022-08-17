#table=read.table("probabilidade de transição 0.7 classificador threshold 0.5 N 1000 p 1",header = T,sep = ";",quote = "\"")

heatmap_matrix_sum=heatmap_array[,,1]
for (k in 2:N_repl){
  heatmap_matrix_sum= heatmap_matrix_sum+heatmap_array[,,k]
}
heatmap_matrix_mean=heatmap_matrix_sum/N_repl


table=as.matrix(heatmap_matrix_mean)
#table=table[,-1]
table1=matrix(NA,11,11)
for(i in 1:11){
  for (j in 1:11){
    table[i,j]=sub(",",".",table[i,j])
    table1[i,j]=as.numeric(as.character(table[i,j]))
  }
}
table1=t(table1)
#levelplot(table1, at=c(0,95,100))
#heatmap(table1)
tab=data.frame()
n=1
l=0
rm(sigma1,sigma2)
sigma1=sigma2=vector()
for (i in rep(seq(1,11,1),11)){
  l=l+1
  sigma1[l]=i
  sigma2[l]=seq(1,11,1)[n]
  if((l%%11)==0)n=n+1
}
resultados=vector()
resultados=table1[,1]
for (k in 2:11){
  resultados[(length(resultados)+1):(length(resultados)+11)]=table1[,k]
}
tab =  as.data.frame(sigma2)
tab[,2] = as.data.frame(seq(from=1,to=2,length=11))
tab[,3] = as.data.frame(table1)
library(ggplot2)
ggplot(tab, aes(sigma1,sigma2, z= resultados)) + geom_tile(aes(fill = resultados)) + 
  theme_bw() + 
  scale_fill_gradient(low="white", high="blue") 

