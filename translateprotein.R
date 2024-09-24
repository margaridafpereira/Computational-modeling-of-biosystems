translateprotein <- function(maiorOrf){
  
  proteina<- translate(s2c(maiorOrf)) #Transforma a string do mRNA em caracteres, função translate traduz o mRNA em aminoácidos 
  
  return(proteina) #Retorna a proteína (aminoácidos)
  
}
proteina<-translateprotein(maiorOrf) #Guarda na variável proteina a tradução do mRNA

