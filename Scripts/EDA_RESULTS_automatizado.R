###################################
#   EDA RESULTS AUTOAMTIZADO      #
###################################
library(readr)
library(dplyr)

#FUNCIONES
add_HKGstats <- function(d){
  n_samples = dim(d)[2]
  d$median = apply(d, 1, median , na.rm = TRUE)
  d$mean = apply(d[,1:n_samples], 1, mean, na.rm = TRUE)
  d$sd = apply(d[,1:n_samples], 1, sd, na.rm = TRUE) # 1-filas
  
  return(d)
}
add_stats <- function(d){
  n_samples = dim(d)[2]
  d$median = apply(d, 1, median , na.rm = TRUE)
  d$mean = apply(d[,1:n_samples], 1, mean, na.rm = TRUE)
  d$sd = apply(d[,1:n_samples], 1, sd, na.rm = TRUE) # 1-filas
  
  
  d = d[order(d$median),]
  d$ranking <- seq.int(nrow(d))
  
  return(d)
}

get_stats <- function(d2){
  #library(dplyr)
  stats_df =  d2 %>% select(ranking, median, mean, sd)
  return(stats_df)
}


# PROGRAMA PRINCIPAL - loop
SEXOS = c("ALL","MALE","FEMALE")
stats = c("CV","IQRmedian","MADmedian")
for (s in 1:length(SEXOS)){
  SEX = SEXOS[s]
  print(SEX)
  for (i in 1:length(stats)){
    stat = stats[i]
    print(stat)
    ### LECTURA FICHEROS:
    print("Lectura de ficheros.")
    if (SEX == "ALL"){
      GPL570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL570_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6244 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6244_",stat,".csv"), header = TRUE, row.names = 1)
      GPL10558 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL10558_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6947 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6947_",stat,".csv"), header = TRUE, row.names = 1)
    }
    if (SEX != "ALL"){
      GPL570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/",SEX,"_GPL570_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6244 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/",SEX,"_GPL6244_",stat,".csv"), header = TRUE, row.names = 1)
      GPL10558 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/",SEX,"_GPL10558_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6947 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/",SEX,"_GPL6947_",stat,".csv"), header = TRUE, row.names = 1)
    }
    ## PARTE 1: RANKING HOUSEKEEPING GENES
    #### PARA LOS VALORES DE LOS HKGS EN TODOS LOS ESTUDIOS:
    
    ## Extraemos los valores para todos los genes (incluimos excepciones del RNA18S):
    #dir.create(paste0("./RESULTS_",stat,"/HKG"))
    #dir.create(paste0("./RESULTS_",stat,"/",SEX))
    print("Creando directorio de resultados HKG.")
    dir.create(paste0("./RESULTS_",stat,"/",SEX,"/HKG"))
    
    print("Interrogando plataformas...")
    GPLs = c("GPL570", "GPL6244", "GPL10558", "GPL6947")
    HKGs = c("HPRT1", "GAPDH","PPIA", "UBC", "RPL19", "RNA18S5")
    for (i in 1:length(HKGs)){
      HKG = HKGs[i]
      if (HKG != "RNA18S5"){
        for (j in 1:length(GPLs)){
          GPLID = GPLs[j]
          df = get(GPLID)  
          which(rownames(df) == HKG)
          df_aux = df[HKG,]
          if (j == 1){HKG_df = df_aux}
          else {HKG_df = cbind(HKG_df, df_aux)}
        }
        #write.csv(HKG_df,paste0("./RESULTS_",stat,"/HKG/",HKG,".csv"), row.names = TRUE)
        ## Creamos el df de los 5 genes con los valores de los 53 estudios:
        if (i == 1) {all_HKG_df = HKG_df} 
        else{all_HKG_df = rbind(all_HKG_df,HKG_df)}
        write.csv(all_HKG_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/5HKG_",stat,".csv"), row.names = TRUE)
      }
      if (HKG == "RNA18S5"){
        GPLs_18S5 = c( "GPL6244", "GPL10558")
        for (r in 1:length(GPLs_18S5)){
          GPLID = GPLs_18S5[r]
          df = get(GPLID)  
          which(rownames(df) == HKG)
          df_aux = df[HKG,]
          if (r == 1){HKG_df = df_aux}
          else {HKG_df = cbind(HKG_df, df_aux)}
        }
        write.csv(HKG_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/RNA18S5_",stat,".csv"), row.names = TRUE)
        RNA18S5 = HKG_df
        }
    }
    ## Calculamos los valores mediano, medio y sd 
    print("Calculando media, mediana y sd.")
    all_HKG_stats = add_HKGstats(all_HKG_df)
    all_HKG_medCV =  all_HKG_stats %>% select(median, mean, sd)
    
    RNA18S_stats = add_HKGstats(RNA18S5)
    RNA18 = RNA18S_stats %>% select(median, mean, sd)
    
    ## Unimos los datos, ordenamos por la mediana y generamos el ranking:
    print("Generando ranking")
    HKG_df = rbind(all_HKG_medCV, RNA18)
    HKG_df = HKG_df[order(HKG_df$median),]
    HKG_df$ranking <- seq.int(nrow(HKG_df))
    
    HKG_df = HKG_df[,c(4,1,2,3)]
    # Guardamos el resultado: 
    print("Escribiendo ranking genes HK.")
    write.csv(HKG_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/HKG_ranking",stat,".csv"), row.names = TRUE)
    
    ## Limpiamos el escritorio de trabajo:
    rm(all_HKG_df, all_HKG_stats, all_HKG_medCV, RNA18, RNA18S, RNA18S_stats,df_aux, HKG_df, GPLID, HKG)
    
    ## PARTE 2: RANKING GLOBAL POR GPL
    ## Creamos directorio para guardar los resultados:
    print("Creando directorio de resultados globalRank")
    dir.create(paste0("./RESULTS_",stat,"/",SEX,"/globalRank"))
    
    ## CALCULAMOS LA SD, MEDIA Y MEDIANA DE LOS CVs
    ## Generamos el ranking basado en la mediana y lo almacenamos en el directorio:
    GPLs = c("GPL570", "GPL6244", "GPL10558", "GPL6947")
    for (i in 1:length(GPLs)){
      GPLID = GPLs[i]
      d = get(GPLID)
      #df_aux = subset (d, select = -c(1))
      df_aux = add_stats(d)
      stats_df = get_stats(df_aux)
      #print(head(stats_df))
      write.csv(stats_df,paste0("./RESULTS_",stat,"/",SEX,"/globalRank/",GPLID,"_rank",stat,".csv"), row.names = TRUE)
      
      #rank_df = get_Medianranking(df_aux)
      #write.csv(rank_df,paste0("./stats/",GPLID,"_rank.csv"), row.names = TRUE)
    }
    ## Limpiamos el escritorio de trabajo:
    rm(d, df_aux,stats_df)
    
    print("Explorando resultados - ¿Dónde caen los genes HK?")
    ## EXPLORAMOS LOS RESULTADOS:
    GPL570_rank <- read.csv(paste0("./RESULTS_",stat,"/",SEX,"/globalRank/GPL570_rank",stat,".csv"), header = TRUE )
    GPL6244_rank <- read.csv(paste0("./RESULTS_",stat,"/",SEX,"/globalRank/GPL6244_rank",stat,".csv"), header = TRUE)
    GPL10558_rank <- read.csv(paste0("./RESULTS_",stat,"/",SEX,"/globalRank/GPL10558_rank",stat,".csv"), header = TRUE)
    GPL6947_rank <- read.csv(paste0("./RESULTS_",stat,"/",SEX,"/globalRank/GPL6947_rank",stat,".csv"), header = TRUE)
    
    ### PARTE 3: ¿En qué posición del ranking global se encuentran los HKGs?
    ### BUSCAMOS LAS POSICIONES DE LOS HKG
    GPLs = c("GPL570", "GPL6244", "GPL10558", "GPL6947")
    nGenes = c(22880, 23307, 31426, 25159)
    #HKGs = c("HPRT1", "GAPDH","PPIA", "UBC", "RPL19")
    for (i in 1:length(GPLs)){
      GPLID = GPLs[i]
      d = get(paste0(GPLID,"_rank"))
      n = nGenes[i]
      text = paste0("\nplatform: ",GPLID," - ",n," genes\n______________________________________________________\n")
      cat(text)
      df_aux = d %>% select(X,ranking, median, mean, sd) %>% filter(X %in% c("HPRT1", "GAPDH","PPIA", "UBC", "RPL19", "RNA18S5"))
      print(df_aux)
      write.csv(df_aux,paste0("./RESULTS_",stat,"/",SEX,"/globalRank/",GPLID,"_positionHKG.csv"))
    }
    }
    
}
