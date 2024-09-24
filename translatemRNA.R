translatemRNA <- function(seqOrf){
  
  seqOrf <- s2c (seqOrf) #Transforma a string do maior ORF em caracteres
  
  mRNA <- seqOrf #mRNA assume a sequência de carateres
  
  for (i in 1:length(seqOrf)){ #Em toda a sequência de caracteres ocorrerão as transcrições apresentadas no respetivo ciclo
    if (seqOrf[i] == 't'){
      mRNA[i]<- 'u' #Nucleótido T passa a U
    }
  }
  mRNA <- c2s (mRNA) #Transforma os caracteres do mRNA em string 
  
  return (mRNA) #Retorna o mRNA
}
mRNA<-translatemRNA(maiorOrf) #Guarda na variável mRNA a transcrição do DNA


