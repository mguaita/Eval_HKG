##############################
#     RANK PRODUCT SCRIPT    #
##############################


## LIBRERIAS
#BiocManager::install("RankProd")
library(RankProd)
library(readr)
library(dplyr)
library(tibble)

## SCRIPT TIENE 2 PARTES:
## Parte 1. Generar los ficheros de entrada
## Parte 2. Función RankProd

## INPUT: 
# GENE | rank_GPL570 | rank_GPL6244 | rank_GPL10558 | rank_GPL6947


## TENEMOS DOS OPCIONES:
# A. utilizar los valores de CV/MAD/IQR - cuanto más pequeño mejor 
# B. utilizar los rankings generados - cuantás arriba en el ranking mejor 

# Ambos datos estan disponibles en los ficheros de salida generados por el script EDA_RESULTS_automatizado.R


### PRIMERO GENERAMOS LOS FICHEROS DE ENTRADA: 
# Generamos las matrices de entrada de la funcion RP con los datos: median y ranking, 
# en total 18 ficheros: 3 CV/IQR/MAD x 3 ALL/MALE/FEMALE x 2 median/ranking 
SEXOS = c("ALL", "MALE", "FEMALE")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      ### LECTURA FICHEROS: generados por el EDA_RESULTS_automated.R - directorios globalRank
      # GEN_ID | valor del estadistico en los GSEs
      GPL570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL570_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      GPL6244 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6244_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      GPL10558 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL10558_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) # row.names = 1)
      GPL6947 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6947_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      
      ## Extraemos la columna que nos interesa - median
      # Contruimos la matriz de entrada para la función RankProd
      GPLs = c("GPL570", "GPL6244", "GPL10558", "GPL6947")
      for (j in 1:length(GPLs)){
        GPLID = GPLs[j]
        df = get(GPLID) 
        df_aux = df[,c("X",option)]
        colnames(df_aux)[2] = paste0(option,GPLID)
        if (j == 1) {RP_df = df_aux}
        if (j > 1){
          # if (GPLID =="GPL6244") { # RNA18S = "RNA18SN5"
          #   # En la plataforma anterior hemos recogido el gen RNA18S como RNA18S5
          #   # Cambiamos el nombre del gen RNA18S en GPL6244 para poder combinarlos en el join
          #   # "RNA18SN5"(GPL6244) a "RNA18S" (GPL570)
          #   df_aux$X[df_aux$X == "RNA18S5"] <- "RNA18S"
          #   #df_aux[29782,1] <- "RNA18S5"
          # }
          # if (GPLID == "GPL10558") { # RNA18S = "RNA18SN5"
          #   # En la plataforma anterior hemos recogido el gen RNA18S como RNA18SN5
          #   # Cambiamos el nombre del gen RNA18S en GPL10558 para poder combinarlos en el join
          #   # "RNA18SN5"(GPL10558) a "RNA18S" (GPL570)
          #   df_aux$X[df_aux$X == "RNA18SN5"] <- "RNA18S"
          #   #df_aux[29782,1] <- "RNA18S5"
          # }
          # Unimos los dfs por la columna X (SYMBOL_ID), es un doble join que introduce NAs a ambos lados, 
          # tendremos genes presentes y ausentes en ambos dfs - introducimos NAs para los genes ausentes
          join_df <- merge(RP_df, df_aux, by= "X", all= TRUE) 
          RP_df = join_df
        }
      }
      # Al finalizar deberiamos tener un df con 41976 genes (41974 si hombres)
      #dir.create(paste0("inputRP_", option))
      write.csv(RP_df,paste0("./RankProd/inputRP_",option,"/",SEX,"_",option,"_inputRP_",stat,".csv"), row.names = FALSE)
    }
  }  
}

## SEGUNDA PARTE: funcion RankProd
setwd("~/Escritorio/Resultados21/RESULTS/RankProd")

SEXOS = c("ALL", "MALE", "FEMALE")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      input_name = paste0("inputRP_median/",SEX,"_median_inputRP_",stat,".csv")
      # Valor mediano del CV de cada gen en cada plataforma
      input_data <- read.csv(input_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      # Posición de los genes en el ranking de cada plataforma:
      #input_data <- read.csv("inputRP_ranking/ALL_ranking_inputRP_CV.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      
      
      ## Definimos nuestros datos como de una sola clase ("1 condición experimental"),
      ## Tenemos que definir una o dos clases.
      ## La función RankProd toma nuestros datos como valores FC
      class = c(1,1,1,1)
      # Extraemos el listado de genes:
      gnames = rownames(input_data)
      # Función RankProd: analisis pair-wise
      # input data
      # class - vector condiciones experimentales, en nuestro caso 1
      # logged - FALSE, nuestros datos no estan en escala logaritmica
      # na.rm - ignorar NAs
      # gene.names - nombre de los genes analizados
      # rand - semilla
      RP_out = RankProducts(input_data, cl = class, logged = FALSE, na.rm =TRUE, gene.names = gnames, plot = FALSE)
      # Mostramos el resultado:
      # num.gene - nº de genes que muestra el resultado, queremos los valores de RP/RSum de todos los genes
      rs = topGene(RP_out, logged = FALSE, num.gene = length(gnames), gene.names = gnames)
      upreg = rs[[1]]
      downreg = rs[[2]]
      
      df = as.data.frame(upreg)
      df_aux = df %>% rownames_to_column("Gene")
      df_aux = df_aux %>% rename("RP" = "RP/Rsum")
      df_aux = df_aux[order(df_aux$RP),]
      df_aux$metaRanking <- seq.int(nrow(df_aux))
      df_aux = df_aux[, c(1,7,3)]
      file_name = paste0("RP_",SEX,"_median",stat,".csv")
      write.csv(df_aux,paste0("./inputRP_median/RESULTS/",file_name), row.names=FALSE)
      
      ### EXPLORACION DE LOS RESULTADOS
      df_aux2 = df_aux %>% select (Gene, metaRanking, RP) %>% filter(Gene %in% c("HPRT1", "GAPDH","PPIA", "UBC", "RPL19", "RNA18S5"))
      print(df_aux2)
      write.csv(df_aux2,paste0("./inputRP_median/RESULTS/RP_",SEX,"_median",stat,"_positionHKG.csv"), row.names = FALSE)
    }
  }
}
### EXPLICACION OUTPUT:
#Two tables of identified genes with:
# gene.index: index of gene in the original data set
# RP/Rsum: Computed rank product/sum for each gene 
# FC:(class1/class2): Expression Fold change of class 1/ class 2. 
# pfp: estimated pfp for each gene if the gene is used as cutoff point 
# P.value: estimated p-value for each gene 
# Table 1 list genes that are up-regulated under class 2, 
# Table 1 ist genes that are down-regulated under class 2
topGene(RP_out, logged = TRUE, num.gene = 10, gene.names = gnames)




### TERCERA PARTE: IDENTIFICAR LOS GENES MÁS ESTABLES
setwd("~/Documentos/TFM/CLUSTER/downloads/RP/RankProd")

SEXOS = c("ALL", "MALE", "FEMALE")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      input_name = paste0("RP_",SEX,"_median",stat,".csv")
      # Valor mediano del CV de cada gen en cada plataforma
      input_data <- read.csv(input_name, header = TRUE)
      if (SEX == "FEMALE"){
        df_aux = input_data[1:1000,]
      }else{
        df_aux = input_data[1:200,]
      }
      write.csv(df_aux,paste0("candidate_genes/RP_",SEX,"_median",stat,"_candidates.csv"), row.names = FALSE)
      
    }
  }
}





#### MUS MUSCULUS ####

SEXOS = c("ALL")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      ### LECTURA FICHEROS: generados por el EDA_RESULTS_automated.R - directorios globalRank
      # GEN_ID | valor del estadistico en los GSEs
      GPL1261 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL1261_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      GPL6246 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6246_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      GPL6885 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6885_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) # row.names = 1)
      GPL6887 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6887_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      GPL16570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL16570_rank",stat,".csv"), header = TRUE, stringsAsFactors = FALSE) #row.names = 1)
      
      ## Extraemos la columna que nos interesa - median
      # Contruimos la matriz de entrada para la función RankProd
      GPLs = c("GPL1261", "GPL6246", "GPL6887", "GPL6885", "GPL16570")
      for (j in 1:length(GPLs)){
        GPLID = GPLs[j]
        df = get(GPLID) 
        df_aux = df[,c("X",option)]
        colnames(df_aux)[2] = paste0(option,GPLID)
        if (j == 1) {RP_df = df_aux}
        if (j > 1){
          # Unimos los dfs por la columna X (SYMBOL_ID), es un doble join que introduce NAs a ambos lados, 
          # tendremos genes presentes y ausentes en ambos dfs - introducimos NAs para los genes ausentes
          join_df <- merge(RP_df, df_aux, by= "X", all= TRUE) 
          RP_df = join_df
        }
      }
      # Al finalizar deberiamos tener un df con 41976 genes (41974 si hombres)
      dir.create(paste0("inputRP_", option))
      write.csv(RP_df,paste0("./RankProd/inputRP_",option,"/",SEX,"_",option,"_inputRP_",stat,".csv"), row.names = FALSE)
    }
  }  
}

## SEGUNDA PARTE: funcion RankProd

SEXOS = c("ALL")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      input_name = paste0("inputRP_median/",SEX,"_median_inputRP_",stat,".csv")
      # Valor mediano del CV de cada gen en cada plataforma
      input_data <- read.csv(input_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      # Posición de los genes en el ranking de cada plataforma:
      #input_data <- read.csv("inputRP_ranking/ALL_ranking_inputRP_CV.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
      
      
      ## Definimos nuestros datos como de una sola clase ("1 condición experimental"),
      ## Tenemos que definir una o dos clases.
      ## La función RankProd toma nuestros datos como valores FC
      class = c(1,1,1,1,1)
      # Extraemos el listado de genes:
      gnames = rownames(input_data)
      # Función RankProd: analisis pair-wise
      # input data
      # class - vector condiciones experimentales, en nuestro caso 1
      # logged - FALSE, nuestros datos no estan en escala logaritmica
      # na.rm - ignorar NAs
      # gene.names - nombre de los genes analizados
      # rand - semilla
      RP_out = RankProducts(input_data, cl = class, logged = FALSE, na.rm =TRUE, gene.names = gnames, plot = FALSE)
      # Mostramos el resultado:
      # num.gene - nº de genes que muestra el resultado, queremos los valores de RP/RSum de todos los genes
      rs = topGene(RP_out, logged = FALSE, num.gene = length(gnames), gene.names = gnames)
      upreg = rs[[1]]
      downreg = rs[[2]]
      
      df = as.data.frame(upreg)
      df_aux = df %>% rownames_to_column("Gene")
      df_aux = df_aux %>% rename("RP" = "RP/Rsum")
      df_aux = df_aux[order(df_aux$RP),]
      df_aux$metaRanking <- seq.int(nrow(df_aux))
      df_aux = df_aux[, c(1,7,3)]
      file_name = paste0("RP_",SEX,"_median",stat,".csv")
      write.csv(df_aux,paste0("inputRP_median/RESULTS/",file_name), row.names=TRUE)
      
      ### EXPLORACION DE LOS RESULTADOS
      df_aux2 = df_aux %>% select(Gene, metaRanking, RP) %>% filter(Gene %in% c("Hprt", "Gapdh","Ppia", "Ubc", "Rpl19", "Rps18"))
      print(df_aux2)
      write.csv(df_aux2,paste0("inputRP_median/RESULTS/RP_",SEX,"_median",stat,"_positionHKG.csv"),row.names = FALSE)
    }
  }
}


SEXOS = c("ALL")
stats = c("CV","IQRmedian","MADmedian")
options = c("median")

for (o in 1:length(options)){
  option = options[o]
  for (s in 1:length(SEXOS)){
    SEX = SEXOS[s]
    for (i in 1:length(stats)){
      stat = stats[i]
      input_name = paste0("RP_",SEX,"_median",stat,".csv")
      # Valor mediano del CV de cada gen en cada plataforma
      input_data <- read.csv(input_name, header = TRUE)
      df_aux = input_data[1:200,]
      write.csv(df_aux,paste0("candidate_genes/RP_",SEX,"_median",stat,"_candidates.csv"), row.names = FALSE)
      
    }
  }
}
