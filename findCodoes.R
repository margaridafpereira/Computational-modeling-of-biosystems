findCodoes<-function(inputseqString){
  
  codoes<-c("atg","taa","tag","tga")
  for (i in 1:4){
    codao<-codoes[i]
    ocorrencias<- matchPattern(codao,inputseqString) #Encontra todas as ocorrências de um determinado codão na sequência
    
    posicoesCodoes<-start(ocorrencias) #Encontra as posições iniciais de todas as ocorrências de cada codão na sequência 
    
    numeroOcorr <- length(posicoesCodoes) #Determina o número de vezes que um certo codão ocorre na sequência
    
    if(i==1){
      posicoes<-posicoesCodoes
      tipos<-rep(codao,numeroOcorr) #Copia para o vetor um determinado codão ("codao"),um determinado número de ocorrências ("numeroOcorr")
      
    }else{
      posicoes<-append(posicoes,posicoesCodoes,after = length(posicoes)) 
      tipos <- append(tipos,rep(codao,numeroOcorr),after = length(tipos))
      #Função "append" adiciona a um determinado vetor,um determinado objeto,numa determinada posição 
    }
  }
  indices <- order(posicoes)
  posicoes<-posicoes[indices]
  tipos<-tipos[indices]
  #Ordena os vetores "posicoes" e "tipos" pela ordem de posições em que ocorrem ao longo da sequência
  
  listaPing<-list(posicoes,tipos) #Cria uma lista que inclui as posições de cada codão na sequência e o respetivo codão
  return(listaPing) 
}



