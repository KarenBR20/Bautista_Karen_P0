##
## 1.secuencia correspondiente de aminoacidos
library(Biostrings)#cargar libreria
x<-readRNAStringSet("C:/Users/EQUIPO/Documents/R/RClasses/ClaseR/Bioinformatica/first.fasta") #subir el archivo
x
translate(x) #traducir el archivo

# 2.Escribe un programa que resuelva los dos problemas que seleccionaste, en ambos, casos, debes buscar una soluci?n sin usar librer?as especializadas y otra con librer?as especializadas.

## 2.1 Counting DNA Nucleotides 
dna <- "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC" #vector donde se guarda la secuencia

#Sin libreria
countdna<- function(){ 
  alfa <-lengths(regmatches(dna, gregexpr("A", dna))) #cantidad de A
  beta <-lengths(regmatches(dna, gregexpr("T", dna)))#cantidad de T
  delta <-lengths(regmatches(dna, gregexpr("C", dna)))#cantidad de C
  omega <-lengths(regmatches(dna, gregexpr("G", dna)))#cantidad de G
  
  return(print(paste(alfa,beta,delta,omega))) #Se imprimen las cantidades de A, T, C y G respectivamente
}

#Con libreria
dna<-DNAStringSet(dna) #el formato del vector tiene que ser DNA strings
letterFrequency(dna, c("A","T","C","G")) #te da las cantidades de las bases nitrogenadas



## 2.2 Computing GC content
raspberry<-readDNAStringSet("https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/raspberry_bushy_dwarf_virus_uid14791/NC_003740.fna") #cargar secuencias
strawberry<-readDNAStringSet("https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/strawberry_latent_ringspot_virus_satellite_rna_uid15155/NC_003848.fna")

#Sin libreria
gc<- function(){ 
  uno <-lengths(regmatches(raspberry, gregexpr("G", raspberry)))#cantidad de G en raspberry
  dos <-lengths(regmatches(raspberry, gregexpr("C", raspberry)))#cantidad de C en raspberry
  tres <-lengths(regmatches(strawberry, gregexpr("G", strawberry)))#cantidad de G en strawberry
  cuatro <-lengths(regmatches(strawberry, gregexpr("C", strawberry)))#cantidad de C en strawberry
  total1<-nchar(raspberry)#numero total de bases en raspberry
  total2<-nchar(strawberry)#numero total de bases en strawberry
  res1 <-(uno+dos)/total1 #frecuencia de GC para raspberry
  res2 <-(tres+cuatro)/total2 #frecuencia de GC para strawberry
  
  if(res1>res2){
    return(print(paste("Raspberry =", res1))) #si raspberry es mas alto se mostrara su frecuencia
  }else if(res2>res1){
    return(print(paste("Strawberry =", res2)))#si strawberry es mas alto se mostrara su frecuencia
  }else{
    return(print(paste("Tienen la misma frecuencia de GC", res1)))#si las frecuencias son iguales mostrara este mensaje
  }
}

#con libreria
porcentaje<- function(){
  one <- letterFrequency(raspberry, c("CG"), as.prob = T)#frecuencia de GC para raspberry
  two <- letterFrequency(strawberry, c("CG"), as.prob = T)#frecuencia de GC para strawberry
  if(one==two){
    return(print(paste("Las frecuencias son iguales"))) #si las frecuencias son iguales mostrara este mensaje
  }else if(one>two){ 
    return(print(paste("Raspberry =", one))) #si raspberry es mas alto se mostrara su frecuencia
  }else{
    return(print(paste("Strawberry =", two))) #si strawberry es mas alto se mostrara su frecuencia
  } 
}