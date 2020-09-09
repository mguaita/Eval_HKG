#LIBRERIAS
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
SEX = "ALL"
stats = c("CV","IQRmedian","MADmedian")
for (i in 1:length(stats)){
  stat = stats[i]
  ### LECTURA FICHEROS:
  GPL1261 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL1261_",stat,".csv"), header = TRUE, row.names = 1)
  GPL6246 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6246_",stat,".csv"), header = TRUE, row.names = 1)
  GPL6887 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6887_",stat,".csv"), header = TRUE, row.names = 1)
  GPL6885 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6885_",stat,".csv"), header = TRUE, row.names = 1)
  GPL16570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL16570_",stat,".csv"), header = TRUE, row.names = 1)
  

  ## PARTE 1: RANKING HOUSEKEEPING GENES
  #### PARA LOS VALORES DE LOS HKGS EN TODOS LOS ESTUDIOS:
  ## Extraemos los valores para todos los genes (incluimos excepciones del RNA18S y Ubc):
  #dir.create(paste0("./RESULTS_",stat,"/HKG"))
 # dir.create(paste0("./RESULTS_",stat,"/",SEX))
  dir.create(paste0("./RESULTS_",stat,"/",SEX,"/HKG"))
  
  GPLs = c("GPL1261", "GPL6246", "GPL6887", "GPL6885", "GPL16570")
  HKGs = c("Hprt", "Gapdh","Ppia", "Ubc", "Rpl19", "Rn18s")
  for (i in 1:length(HKGs)){
    HKG = HKGs[i]
    if (HKG != "Rn18s" && HKG != "Ubc"){
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
    if (HKG == "Rn18s"){
      RNA18s_GPLs = c("GPL1261", "GPL6246", "GPL6887")
      for (j in 1:length(RNA18s_GPLs)){
        GPLID = RNA18s_GPLs[j]
        df = get(GPLID)  
        which(rownames(df) == HKG)
        df_aux = df[HKG,]
        if (j == 1){Rn18s_df = df_aux}
        else {Rn18s_df = cbind(Rn18s_df, df_aux)}
      }
      write.csv(Rn18s_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/Rn18s_",stat,".csv"), row.names = TRUE)
    }
    if (HKG == "Ubc"){
      Ubc_GPLs = c("GPL1261", "GPL6246", "GPL6885", "GPL16570")
      for (u in 1:length(Ubc_GPLs)){
          GPLID =  Ubc_GPLs[u]
          df = get(GPLID)  
          which(rownames(df) == HKG)
          df_aux = df[HKG,]
          if (u == 1){Ubc_df = df_aux}
          else {Ubc_df = cbind(Ubc_df, df_aux)}
      }
      write.csv(Ubc_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/Ubc_",stat,".csv"), row.names = TRUE)
      
      }
  }
  ## Calculamos los valores mediano, medio y sd 
  all_HKG_stats = add_HKGstats(all_HKG_df)
  all_HKG_medCV =  all_HKG_stats %>% select(median, mean, sd)
  
  RNA18S_stats = add_HKGstats(Rn18s_df)
  RNA18 = RNA18S_stats %>% select(median, mean, sd)
  
  Ubc_stats = add_HKGstats(Ubc_df)
  UBC = Ubc_stats %>% select(median, mean, sd)
  
  ## Unimos los datos, ordenamos por la mediana y generamos el ranking:
  HKG_df = rbind(all_HKG_medCV, UBC, RNA18)
  HKG_df = HKG_df[order(HKG_df$median),]
  HKG_df$ranking <- seq.int(nrow(HKG_df))
  
  HKG_df = HKG_df[,c(4,1,2,3)]
  # Guardamos el resultado: 
  write.csv(HKG_df,paste0("./RESULTS_",stat,"/",SEX,"/HKG/HKG_ranking",stat,".csv"), row.names = TRUE)
  
  ## Limpiamos el escritorio de trabajo:
  rm(all_HKG_df, all_HKG_stats, all_HKG_medCV, RNA18, RNA18S, RNA18S_stats, Ubc_stats, UBC, df_aux, HKG_df, GPLID, HKG)
  
  ## PARTE 2: RANKING GLOBAL POR GPL
  ## Creamos directorio para guardar los resultados:
  dir.create(paste0("./RESULTS_",stat,"/",SEX,"/globalRank"))
  
  ## CALCULAMOS LA SD, MEDIA Y MEDIANA DE LOS CVs
  ## Generamos el ranking basado en la mediana y lo almacenamos en el directorio:
  GPLs = c("GPL1261", "GPL6246", "GPL6887", "GPL6885", "GPL16570")
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
  
  ## EXPLORAMOS LOS RESULTADOS:
  
  GPL1261_rank  <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL1261_rank",stat,".csv"), header = TRUE )
  GPL6246_rank  <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6246_rank",stat,".csv"), header = TRUE )
  GPL6887_rank  <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6887_rank",stat,".csv"), header = TRUE )
  GPL6885_rank  <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL6885_rank",stat,".csv"), header = TRUE )
  GPL16570_rank  <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/globalRank/GPL16570_rank",stat,".csv"), header = TRUE )
  

  ### PARTE 3: ¿En qué posición del ranking global se encuentran los HKGs?
  ### BUSCAMOS LAS POSICIONES DE LOS HKG
  GPLs = c("GPL1261", "GPL6246", "GPL6887", "GPL6885", "GPL16570")
  nGenes = c(21496,24213,30866,18120,24647)
  #HKGs = c("HPRT1", "GAPDH","PPIA", "UBC", "RPL19")
  for (i in 1:length(GPLs)){
    GPLID = GPLs[i]
    d = get(paste0(GPLID,"_rank"))
    n = nGenes[i]
    text = paste0("\nplatform: ",GPLID," - ",n," genes\n______________________________________________________\n")
    cat(text)
    df_aux = d %>% select(X,ranking, median, mean, sd) %>% filter(X %in% c("Hprt", "Gapdh","Ppia", "Ubc", "Rpl19", "Rn18s"))
    print(df_aux)
    write.csv(df_aux,paste0("./RESULTS_",stat,"/",SEX,"/globalRank/",GPLID,"_positionHKG.csv"))
  }
}
