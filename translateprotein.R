translateprotein <- function(maiorOrf){
  
  proteina<- translate(s2c(maiorOrf)) #Transforma a string do mRNA em caracteres, fun��o translate traduz o mRNA em amino�cidos 
  
  return(proteina) #Retorna a prote�na (amino�cidos)
  
}
proteina<-translateprotein(maiorOrf) #Guarda na vari�vel proteina a tradu��o do mRNA

