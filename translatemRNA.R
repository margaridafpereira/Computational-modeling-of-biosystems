translatemRNA <- function(seqOrf){
  
  seqOrf <- s2c (seqOrf) #Transforma a string do maior ORF em caracteres
  
  mRNA <- seqOrf #mRNA assume a sequ�ncia de carateres
  
  for (i in 1:length(seqOrf)){ #Em toda a sequ�ncia de caracteres ocorrer�o as transcri��es apresentadas no respetivo ciclo
    if (seqOrf[i] == 't'){
      mRNA[i]<- 'u' #Nucle�tido T passa a U
    }
  }
  mRNA <- c2s (mRNA) #Transforma os caracteres do mRNA em string 
  
  return (mRNA) #Retorna o mRNA
}
mRNA<-translatemRNA(maiorOrf) #Guarda na vari�vel mRNA a transcri��o do DNA


