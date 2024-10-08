slidingwindowGC<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  
  n<-length(starts) #Determina o comprimento do vetor "starts"
  
  chunkGCs<-numeric(n) #Cria um vetor com o mesmo comprimento que o vetor "starts" mas, apenas contendo zeros 
  
  for(i in 1:n){
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequ�ncia para ser analisada
    
    chunkGC<-GC(chunk) #Na sequ�ncia anterior conta a frequ�ncia dos nucle�tidos "C" e "G" e  soma-os
    
    chunkGCs[i]<-chunkGC #Guarda em posi��es sucessivas o valor da frequ�ncia do conte�do "GC"
    
  }
  plot(starts,chunkGCs,type="l",xlab="Posi��o inicial do nucle�tido",ylab="Conte�do GC") #Legenda os eixos x e y
  
  conteudoGC<-GC(pinguimseq) #Conte�do "GC" na sequ�ncia
  
  print(conteudoGC) #Faz print do conte�do "GC"
}

