findOrfs<-function(inputseq){
  inputseq<-c2s(inputseq) #Transforma a sequência de caracteres numa string
  
  pinguimLista<-findCodoes(inputseq) #Vai buscar informação à função findCodoes
  
  posicoes<- pinguimLista[[1]] #Vai buscar as posições à função findCodoes
  
  codoes<-pinguimLista[[2]] #Vai buscar os codões à função findCodoes
  
  posIniciais <- numeric() #Vetor de números
  posFinais<- numeric() #Vetor de números
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
        #Como a cond1 e a cond2 passaram a 1, encontrado codão de iniciação
      }
    }else if(cond1==1 && cond2==1){
      if(codoes[i]=="taa"||codoes[i]=="tga"||codoes[i]=="tag"){
        if((((posicoes[i])-posIniciais[length(posIniciais)])%%3)==0){
          posFinais<-append(posFinais,posicoes[i]+2)
          cond1<-0
          cond2<-0
          #A cond1 e a cond2 passam a 0 se for encontrado um codão de finalização
        }else{
          cond1<- 1
          cond2<- 1
          #A cond1 e a cond2 continuam a 1 se não for encontrado um codão de finalização
        }
      }
    }
  }
  for(i in 1:length(posFinais)){
    a<-substring(inputseq,posIniciais[i],posFinais[i])
    subSeq<-append(subSeq,a)
    #Quebra a sequência de DNA em ORFs
  }
  for(i in 1:length(subSeq)){
    if(nchar(subSeq[i])>maior){
      maior<-nchar(subSeq[i])
      posMaior<-i
      #Na lista de ORF obtida procura qual o maior
    }
  }
  maiorOrf<-subSeq[posMaior] #Guarda a sequência do maior ORF
  return(maiorOrf)
}
maiorOrf<-findOrfs(pinguimseq) #Guarda na variável maiorOrf a sequência do maior ORF




