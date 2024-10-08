slidingwindowAT<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  
  n<-length(starts) #Determina o comprimento do vetor "starts" mas, apenas contendo zeros
  
  chunkATs<-numeric(n) #Cria um vetor com o mesmo comprimento que o vetor "starts" mas, apenas contendo zeros
  
  for(i in 1:n)
  {
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequ�ncia para ser analisada
    
    chunkAT<-count(chunk,1,freq = T)["a"]+count(chunk,1,freq = T)["t"] #Na sequ�ncia anterior conta a frequ�ncia dos nucle�tidos "A" e "T" e soma-os
    
    chunkATs[i]<-chunkAT #Guarda em posi��es sucessivas o valor da frequ�ncia do conte�do "AT"
  }
  plot(starts,chunkATs,type="l",xlab="Posi��o inicial do nucle�tido",ylab="Conte�do AT",col="black") # Atribui��o da cor � linha do gr�fico e legenda dos eixos
}

