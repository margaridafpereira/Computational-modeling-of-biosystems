slidingwindowAT<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  
  n<-length(starts) #Determina o comprimento do vetor "starts" mas, apenas contendo zeros
  
  chunkATs<-numeric(n) #Cria um vetor com o mesmo comprimento que o vetor "starts" mas, apenas contendo zeros
  
  for(i in 1:n)
  {
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequência para ser analisada
    
    chunkAT<-count(chunk,1,freq = T)["a"]+count(chunk,1,freq = T)["t"] #Na sequência anterior conta a frequência dos nucleótidos "A" e "T" e soma-os
    
    chunkATs[i]<-chunkAT #Guarda em posições sucessivas o valor da frequência do conteúdo "AT"
  }
  plot(starts,chunkATs,type="l",xlab="Posição inicial do nucleótido",ylab="Conteúdo AT",col="black") # Atribuição da cor à linha do gráfico e legenda dos eixos
}

