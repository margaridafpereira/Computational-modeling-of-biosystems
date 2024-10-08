findCodoes<-function(inputseqString){
  
  codoes<-c("atg","taa","tag","tga")
  for (i in 1:4){
    codao<-codoes[i]
    ocorrencias<- matchPattern(codao,inputseqString) #Encontra todas as ocorr�ncias de um determinado cod�o na sequ�ncia
    
    posicoesCodoes<-start(ocorrencias) #Encontra as posi��es iniciais de todas as ocorr�ncias de cada cod�o na sequ�ncia 
    
    numeroOcorr <- length(posicoesCodoes) #Determina o n�mero de vezes que um certo cod�o ocorre na sequ�ncia
    
    if(i==1){
      posicoes<-posicoesCodoes
      tipos<-rep(codao,numeroOcorr) #Copia para o vetor um determinado cod�o ("codao"),um determinado n�mero de ocorr�ncias ("numeroOcorr")
      
    }else{
      posicoes<-append(posicoes,posicoesCodoes,after = length(posicoes)) 
      tipos <- append(tipos,rep(codao,numeroOcorr),after = length(tipos))
      #Fun��o "append" adiciona a um determinado vetor,um determinado objeto,numa determinada posi��o 
    }
  }
  indices <- order(posicoes)
  posicoes<-posicoes[indices]
  tipos<-tipos[indices]
  #Ordena os vetores "posicoes" e "tipos" pela ordem de posi��es em que ocorrem ao longo da sequ�ncia
  
  listaPing<-list(posicoes,tipos) #Cria uma lista que inclui as posi��es de cada cod�o na sequ�ncia e o respetivo cod�o
  return(listaPing) 
}



