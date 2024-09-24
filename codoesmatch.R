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

