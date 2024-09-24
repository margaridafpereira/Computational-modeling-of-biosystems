slidingwindowGC<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  
  n<-length(starts) #Determina o comprimento do vetor "starts"
  
  chunkGCs<-numeric(n) #Cria um vetor com o mesmo comprimento que o vetor "starts" mas, apenas contendo zeros 
  
  for(i in 1:n){
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequência para ser analisada
    
    chunkGC<-GC(chunk) #Na sequência anterior conta a frequência dos nucleótidos "C" e "G" e  soma-os
    
    chunkGCs[i]<-chunkGC #Guarda em posições sucessivas o valor da frequência do conteúdo "GC"
    
  }
  plot(starts,chunkGCs,type="l",xlab="Posição inicial do nucleótido",ylab="Conteúdo GC") #Legenda os eixos x e y
  
  conteudoGC<-GC(pinguimseq) #Conteúdo "GC" na sequência
  
  print(conteudoGC) #Faz print do conteúdo "GC"
}

