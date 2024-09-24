slidingwindownucleotidos<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  n<-length(starts) #Determina o comprimento do vetor "starts" mas, apenas contendo zeros
  chunkAs<-numeric(n) 
  chunkTs<-numeric(n) 
  chunkGs<-numeric(n) 
  chunkCs<-numeric(n) 
  # Cria vetores com o mesmo comprimento que o vector "starts" mas, apenas contendo zeros
  
  for(i in 1:n){
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequência para ser analisada
    chunkA<-count(chunk,1,freq = T)["a"] 
    chunkT<-count(chunk,1,freq = T)["t"] 
    chunkG<-count(chunk,1,freq = T)["g"] 
    chunkC<-count(chunk,1,freq = T)["c"] 
    # Armazenam a frequência de cada um dos nucleótidos
    
    chunkAs[i]<-chunkA 
    chunkTs[i]<-chunkT 
    chunkGs[i]<-chunkG 
    chunkCs[i]<-chunkC 
    # Guardam em posições sucessivas o valor da frequência de cada um dos nucleótidos
    
  }
  lim2<-max(chunkAs,chunkTs,chunkGs,chunkCs)
  lim1<-min(chunkAs,chunkTs,chunkGs,chunkCs)
  plot(starts,chunkAs,ylim=range(c(lim1,lim2)),type="l",xlab="Posição inicial do nucleótido",ylab="Conteúdo dos nucleótidos",col="blue") 
  par(new=TRUE)
  plot(starts,chunkTs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="red") 
  par(new=TRUE)
  plot(starts,chunkGs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="yellow") 
  par(new=TRUE)
  plot(starts,chunkCs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="green") 
  legend("bottomleft",legend=c("A","T","G","C"),text.col = c("blue","red","yellow","green"),box.lty = 0,horiz = F)
  
  # Atribuição da respetiva cor a cada uma das  linhas do gráfico e legenda dos eixos
}

