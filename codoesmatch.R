codoesmatch<-function(pinguimseq){
  codoesin<-c("atg") #atg � o cod�o start que vai ser armazenado pelo codoesin
  
  codao<-codoesin[1] #cria o vetor codao para conter os cod�es start
  
  ocorrenciasin<-matchPattern(codao,c2s(pinguimseq)) #pesquisa na sequ�ncia pinguimseq o "atg"
  
  tamanho<-length(ocorrenciasin) #vari�vel tamanho armazena o n�mero de ocorr�ncias do cod�o start
  
  print("N�mero de ocorr�ncias do cod�o de start:")
  print(tamanho) #print do n�mero de cod�es start
  
  codoesfim<-c("taa","tag","tga") #"taa","tag" e "tga" s�o os cod�es stop que ser�o armazenados pelo codoesfim
  tamanhototal<-0
  for(i in 1:3)
  {
    codaofim<-codoesfim[i] #cria um vetor para os cod�es stop
    
    ocorrenciasfim<-matchPattern(codaofim,c2s(pinguimseq)) #pesquisa na sequ�ncia as 3 combina��es
    
    tamanhofinal<-length(ocorrenciasfim) #armazena o n�mero de ocorr�ncias dos cod�es de stop
    
    tamanhototal<-tamanhototal+tamanhofinal #calcula o n�mero de cod�es stop
  }
  print("N�mero de ocorr�ncias de cod�es stop:")
  print(tamanhototal) #print do n�mero de cod�es stop
}

