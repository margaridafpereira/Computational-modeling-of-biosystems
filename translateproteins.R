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


