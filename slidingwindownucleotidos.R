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
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequ�ncia para ser analisada
    chunkA<-count(chunk,1,freq = T)["a"] 
    chunkT<-count(chunk,1,freq = T)["t"] 
    chunkG<-count(chunk,1,freq = T)["g"] 
    chunkC<-count(chunk,1,freq = T)["c"] 
    # Armazenam a frequ�ncia de cada um dos nucle�tidos
    
    chunkAs[i]<-chunkA 
    chunkTs[i]<-chunkT 
    chunkGs[i]<-chunkG 
    chunkCs[i]<-chunkC 
    # Guardam em posi��es sucessivas o valor da frequ�ncia de cada um dos nucle�tidos
    
  }
  lim2<-max(chunkAs,chunkTs,chunkGs,chunkCs)
  lim1<-min(chunkAs,chunkTs,chunkGs,chunkCs)
  plot(starts,chunkAs,ylim=range(c(lim1,lim2)),type="l",xlab="Posi��o inicial do nucle�tido",ylab="Conte�do dos nucle�tidos",col="blue") 
  par(new=TRUE)
  plot(starts,chunkTs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="red") 
  par(new=TRUE)
  plot(starts,chunkGs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="yellow") 
  par(new=TRUE)
  plot(starts,chunkCs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="green") 
  legend("bottomleft",legend=c("A","T","G","C"),text.col = c("blue","red","yellow","green"),box.lty = 0,horiz = F)
  
  # Atribui��o da respetiva cor a cada uma das  linhas do gr�fico e legenda dos eixos
}

