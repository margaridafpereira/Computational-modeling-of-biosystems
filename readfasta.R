pinguim<-read.fasta(file = "NC_021474.1.fasta") #leitura de arquivo no formato FASTA
pinguimseq<-pinguim[[1]] #pinguimseq contém a sequência do genoma - criação de um vetor
pinguimSeqString <- c2s(pinguimseq)