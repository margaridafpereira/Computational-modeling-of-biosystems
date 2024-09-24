#PROGRAMA

#1. Pinguim (Spheniscidae)

#2.1. Pinguim - Código de acesso: NC_021474.1


library(BiocGenerics)
library(BiocManager)
library(BiocVersion)
library(Biostrings)
library(seqinr)

#2.2.


pinguim<-read.fasta(file = "NC_021474.1.fasta") #leitura de arquivo no formato FASTA
pinguimseq<-pinguim[[1]] #pinguimseq contém a sequência do genoma - criação de um vetor
pinguimSeqString <- c2s(pinguimseq)
length(pinguimseq)
write.fasta(name="pinguim", sequences = pinguim, file.out = "pinguim.fasta") #gravação num ficheiro FASTA


#3.1.1.


count(pinguimseq,1)
count(pinguimseq,1,freq = T)
barplot(table(pinguimseq))
plot(table(pinguimseq))
pie(table(pinguimseq))


#3.1.2.


count(pinguimseq,3)
count(pinguimseq,3,freq = T)


#3.1.3.


codoesmatch<-function(pinguimseq){
  codoesin<-c("atg") #atg é o codão start que vai ser armazenado pelo codoesin
  codao<-codoesin[1] #cria o vetor codao para conter os codões start
  ocorrenciasin<-matchPattern(codao,c2s(pinguimseq)) #pesquisa na sequência pinguimseq o "atg"
  tamanho<-length(ocorrenciasin) #variável tamanho armazena o número de ocorrências do codão start
  print("Número de ocorrências do codão de start:")
  print(tamanho) #print do número de codões start
  codoesfim<-c("taa","tag","tga") #"taa","tag" e "tga" são os codões stop que serão armazenados pelo codoesfim
  tamanhototal<-0
  for(i in 1:3)
  {
    codaofim<-codoesfim[i] #cria um vetor para os codões stop
    ocorrenciasfim<-matchPattern(codaofim,c2s(pinguimseq)) #pesquisa na sequência as 3 combinações
    tamanhofinal<-length(ocorrenciasfim) #armazena o número de ocorrências dos codões de stop
    tamanhototal<-tamanhototal+tamanhofinal #calcula o número de codões stop
  }
  print("Número de ocorrências de codões stop:")
  print(tamanhototal) #print do número de codões stop
}


#3.1.4.


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



#3.1.4.1.


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


#3.1.4.2.


translateprotein <- function(maiorOrf){
  proteina<- translate(s2c(maiorOrf)) #Transforma a string do mRNA em caracteres, função translate traduz o mRNA em aminoácidos 
  return(proteina) #Retorna a proteína (aminoácidos)
}
proteina<-translateprotein(maiorOrf) #Guarda na variável proteina a tradução do mRNA


translateproteins <- function(maiorOrf){
  codoes<-sapply(seq(1,nchar(mRNA),3), function(i)substring(mRNA,i,i+2)) #Divide o mRNA em sequências de 3
  n<-length(codoes) #Determina o comprimento do vetor "codoes"
  proteinas<-character(n) #Cria um vetor vazio com o mesmo comprimento que o vetor "codoes"
  
  for(i in 1:length(codoes)){
    
    if(codoes[i] == "ggu" || codoes[i] =="ggc"|| codoes[i] == "gga" || codoes[i] =="ggg"){
      
      proteinas[i]<-"Gli" 
    }
    
    if(codoes[i] == "uuu" || codoes[i] == "uuc"){
      
      proteinas[i]<-"Fen"
    }
    
    if(codoes[i] == "uua" || codoes[i] =="uug"||codoes[i] == "cuu"||codoes[i] =="cuc"||codoes[i] =="cua"||codoes[i] =="cug"){
      
      proteinas[i]<-"Leuc"
    }
    
    if(codoes[i] == "auu" || codoes[i] =="auc"||codoes[i] == "aua"){
      
      proteinas[i]<-"Isol"
    }
    
    if(codoes[i] == "aug" ){
      
      proteinas[i]<-"Met"
    }
    
    if(codoes[i] == "guu" || codoes[i] =="guc"||codoes[i] == "gua"|| codoes[i] == "gug"){
      
      proteinas[i]<-"Val"  
    }
    
    if(codoes[i] == "ucu" || codoes[i] =="ucc"||codoes[i] == "uca"||codoes[i] == "ucg"||codoes[i] == "agu"||codoes[i] == "agc"){
      
      proteinas[i]<-"Ser"  
    }
    
    if(codoes[i] == "ccu" || codoes[i] =="cca"||codoes[i] == "ccc"||codoes[i] == "ccg"){
      
      proteinas[i]<-"Prol"
    } 
    
    if(codoes[i] == "acu" || codoes[i] =="acc"||codoes[i] == "aca" || codoes[i] == "acg"){
      
      proteinas[i]<-"Treo" 
    }
    
    if(codoes[i] == "gcu" || codoes[i] =="gcc"||codoes[i] == "gca"||codoes[i] == "gcg"){
      
      proteinas[i]<-"Ala" 
    }
    
    if(codoes[i] == "cgu" || codoes[i] =="cga"|| codoes[i] == "cgc" || codoes[i] =="cgg"|| codoes[i] == "aga" || codoes[i] =="agg"){
      
      proteinas[i]<-"Arg" 
    } 
    
    if(codoes[i] == "uau" || codoes[i] =="uac"){
      
      proteinas[i]<-"Tir" 
    }
    
    if(codoes[i] == "cau" || codoes[i] =="cac"){
      
      proteinas[i]<-"Hist" 
    }
    
    if(codoes[i] == "caa" || codoes[i] =="cag"){
      
      proteinas[i]<-"Glut" 
    }
    
    if(codoes[i] == "aau" || codoes[i] =="aac"){
      
      proteinas[i]<-"Asp" 
    }
    
    if(codoes[i] == "aaa" || codoes[i] =="aag"){
      
      proteinas[i]<-"Lis" 
    }
    
    if(codoes[i] == "gau" || codoes[i] =="gac"){
      
      proteinas[i]<-"AcAsp" 
    } 
    
    if(codoes[i] == "ugg" ){
      
      proteinas[i]<-"Trip" 
    }
    
    if(codoes[i] == "gaa" || codoes[i] =="gag"){
      
      proteinas[i]<-"AcGlut" 
    }
    
    if(codoes[i] == "ugu" || codoes[i] =="ugc"){
      
      proteinas[i]<-"Cist" 
    }
    
    #Associa a cada codão o respetivo aminoácido
  }
  return(proteinas)
}
proteinas<-translateproteins(maiorOrf) #Guarda na variável proteinas a tradução do mRNA


#3.1.5.1.


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


#3.1.5.2.


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


#3.1.5.3.


slidingwindownucleotidos<-function(windowsize,inputseq)
{
  starts<-seq(1,length(inputseq)-windowsize)
  n<-length(starts) #Determina o comprimento do vetor "starts" mas, apenas contendo zeros
  chunkAs<-numeric(n) 
  chunkTs<-numeric(n) 
  chunkGs<-numeric(n) 
  chunkCs<-numeric(n) 
  # Cria vetores com o mesmo comprimento que o vector "starts" mas, apenas contendo zeros
  for(i in 1:n){
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)] #Guarda parte da sequência para ser analisada
    chunkA<-count(chunk,1,freq = T)["a"] 
    chunkT<-count(chunk,1,freq = T)["t"] 
    chunkG<-count(chunk,1,freq = T)["g"] 
    chunkC<-count(chunk,1,freq = T)["c"] 
    # Armazenam a frequência de cada um dos nucleótidos
    chunkAs[i]<-chunkA 
    chunkTs[i]<-chunkT 
    chunkGs[i]<-chunkG 
    chunkCs[i]<-chunkC 
    # Guardam em posições sucessivas o valor da frequência de cada um dos nucleótidos
  }
  lim2<-max(chunkAs,chunkTs,chunkGs,chunkCs)
  lim1<-min(chunkAs,chunkTs,chunkGs,chunkCs)
  plot(starts,chunkAs,ylim=range(c(lim1,lim2)),type="l",xlab="Posição inicial do nucleótido",ylab="Conteúdo dos nucleótidos",col="blue") 
  par(new=TRUE)
  plot(starts,chunkTs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="red") 
  par(new=TRUE)
  plot(starts,chunkGs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="yellow") 
  par(new=TRUE)
  plot(starts,chunkCs,ylim=range(c(lim1,lim2)),type="l",xlab="",ylab="",col="green") 
  legend("bottomleft",legend=c("A","T","G","C"),text.col = c("blue","red","yellow","green"),box.lty = 0,horiz = F)
  # Atribuição da respetiva cor a cada uma das  linhas do gráfico e legenda dos eixos
}
