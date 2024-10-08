findOrfs<-function(inputseq){
  inputseq<-c2s(inputseq) #Transforma a sequ�ncia de caracteres numa string
  
  pinguimLista<-findCodoes(inputseq) #Vai buscar informa��o � fun��o findCodoes
  
  posicoes<- pinguimLista[[1]] #Vai buscar as posi��es � fun��o findCodoes
  
  codoes<-pinguimLista[[2]] #Vai buscar os cod�es � fun��o findCodoes
  
  posIniciais <- numeric() #Vetor de n�meros
  posFinais<- numeric() #Vetor de n�meros
  subSeq<- character() #Vetor de caracteres
  cond1 <-0
  cond2<-0
  maior<-0
  posMaior<-0
  for (i in 1:length(codoes)){
    if(cond1==0 && cond2==0){
      if(codoes[i]=="atg"){
        posIniciais<-append(posIniciais,posicoes[i])
        cond1<-1
        cond2<-1
        #Como a cond1 e a cond2 passaram a 1, encontrado cod�o de inicia��o
      }
    }else if(cond1==1 && cond2==1){
      if(codoes[i]=="taa"||codoes[i]=="tga"||codoes[i]=="tag"){
        if((((posicoes[i])-posIniciais[length(posIniciais)])%%3)==0){
          posFinais<-append(posFinais,posicoes[i]+2)
          cond1<-0
          cond2<-0
          #A cond1 e a cond2 passam a 0 se for encontrado um cod�o de finaliza��o
        }else{
          cond1<- 1
          cond2<- 1
          #A cond1 e a cond2 continuam a 1 se n�o for encontrado um cod�o de finaliza��o
        }
      }
    }
  }
  for(i in 1:length(posFinais)){
    a<-substring(inputseq,posIniciais[i],posFinais[i])
    subSeq<-append(subSeq,a)
    #Quebra a sequ�ncia de DNA em ORFs
  }
  for(i in 1:length(subSeq)){
    if(nchar(subSeq[i])>maior){
      maior<-nchar(subSeq[i])
      posMaior<-i
      #Na lista de ORF obtida procura qual o maior
    }
  }
  maiorOrf<-subSeq[posMaior] #Guarda a sequ�ncia do maior ORF
  return(maiorOrf)
}
maiorOrf<-findOrfs(pinguimseq) #Guarda na vari�vel maiorOrf a sequ�ncia do maior ORF




