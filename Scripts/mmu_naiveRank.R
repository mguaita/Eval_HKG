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
SEXOS = c("ALL")
stats = c("CV","IQRmedian","MADmedian")

for (s in 1:length(SEXOS)){
  SEX = SEXOS[s]
  for (i in 1:length(stats)){
    stat = stats[i]
    
    ### LECTURA FICHEROS:
    if (SEX == "ALL"){
      GPL1261 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL1261_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6246 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6246_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6887 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6887_",stat,".csv"), header = TRUE, row.names = 1)
      GPL6885 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL6885_",stat,".csv"), header = TRUE, row.names = 1)
      GPL16570 <- read.csv(paste0("RESULTS_",stat,"/",SEX,"/GPL16570_",stat,".csv"), header = TRUE, row.names = 1)
    }
      
    ## PARTE 4: METANALISIS NAIVE - sin asignar pesos
    # Recogemos todos los identificadores SYMBOL de todas las plataformas, e interrogamos todos los genes.
    # Para cada gen que somos capaces de identificar, recogemos los valores del estadistico en todos los estudios 
    # en los que aparezca, al generar el ranking general nos quedamos con el valor mediano.
    
    # GPL570 - 22880 genes
    # GPL6244 - 23307 genes
    # GPL10558 - 31426 genes
    # GPL6947 - 25159 genes
    # Total - 102772 identificadores SYMBOL
    
    # Extra los SYMBOL de cada plataforma y los almacenamos en una lista:
    GPLs = c("GPL1261", "GPL6246", "GPL6887", "GPL6885", "GPL16570")
    SYMBOL_ID = c()
    for (j in 1:length(GPLs)){
      GPLID = GPLs[j]
      df = get(GPLID) 
      SYMBOL = rownames(df)
      SYMBOL_ID = append(SYMBOL_ID, SYMBOL)
    }
    rm(df,j, GPLID,SYMBOL)
    
    # Eliminamos los ID repetidos: 102772 identificadores -> 41976 identificadores
    SYMBOL_ID = unique(SYMBOL_ID) 
    MEDIAN_VALUE = c() # Vector par almacenar los valores medianos
    nGSEs = c() # vector para almacenar el nº de estudios en los que aparece cada gen
    # Para cada estadistico: CV, IQR, MAD
    # Interrogamos cada gen de la lista SYMBOL_ID para extraer sus valores en cada plataforma:
    # 1. Extraemos sus valores de los estudios de cada plataforma
    # 2. Unimos los resultados en un unico vector/fila de todos los estudios en los que aparece (GEN_stat_vector)
    # 3. Obtenemos el valor mediano y anotamos el nº de estudios en los que aparece (length(GEN_stat_vector))
    # 4. Almacenamos el valor mediano en un vector de medianas (MEDIAN_VALUE)
    # 5. Construimos un df con 3 columnas, 1. SYMBOL_ID, 2. MEDIAN_VALUE, 3. nGSE_GEN
    # 6. Generamos un ranking ordenando por este valor mediano (columna3)
    for (i in 1:length(SYMBOL_ID)){
      GEN = SYMBOL_ID[i]
      GEN_stat_vector = c()
      # Sabemos que RNA18S5 tiene dos alias y solo aparece en 2 plataformas
      # Tratamos sus resultados aparte porque tiene 2 "entradas" en lista de genes
      # Se interroga 2 veces
      if (GEN != "Rn18s" && GEN != "Ubc"){ # El resto de genes se interrogan para las 4 GPLs
        for (j in 1:length(GPLs)){
          GPLID = GPLs[j]
          df = get(GPLID)  
          which(rownames(df) == GEN)
          df_aux = df[GEN,]
          # Convertimos el df en un vector
          v_aux = as.numeric(as.vector(df[GEN,]))
          GEN_stat_vector = c(GEN_stat_vector, v_aux)}
        # Extraidos los valores para un gen:
        # anotamos el nº de studios en los que aparece y calculamos el valor mediano
        # Nº de estudios en los que aparece = longitud vector almacenamiento menos los NAs
        nGSE_GEN = length(GEN_stat_vector) - sum(is.na(GEN_stat_vector))
        nGSEs = c(nGSEs,nGSE_GEN)
        # Obtenemos el valor mediano (sin contar NAs) y lo almacenamos:
        MEDIAN = median(GEN_stat_vector, na.rm = TRUE)
        MEDIAN_VALUE = c(MEDIAN_VALUE, MEDIAN)
      }else if(GEN == "Rn18s"){
        RNA18s_GPLs = c("GPL1261", "GPL6246", "GPL6887")
        for (j in 1:length(RNA18s_GPLs)){
          GPLID = RNA18s_GPLs[j]
          df = get(GPLID)  
          which(rownames(df) == GEN)
          df_aux = df[GEN,]
          # Convertimos el df en un vector
          v_aux = as.numeric(as.vector(df[GEN,]))
          GEN_stat_vector = c(GEN_stat_vector, v_aux)}
        # Extraidos los valores para un gen:
        # anotamos el nº de studios en los que aparece y calculamos el valor mediano
        # Nº de estudios en los que aparece = longitud vector almacenamiento menos los NAs
        nGSE_GEN = length(GEN_stat_vector) - sum(is.na(GEN_stat_vector))
        nGSEs = c(nGSEs,nGSE_GEN)
        # Obtenemos el valor mediano (sin contar NAs) y lo almacenamos:
        MEDIAN = median(GEN_stat_vector, na.rm = TRUE)
        MEDIAN_VALUE = c(MEDIAN_VALUE, MEDIAN)
      }else if (GEN == "Ubc"){
        Ubc_GPLs = c("GPL1261", "GPL6246", "GPL6885", "GPL16570")
        for (u in 1:length(Ubc_GPLs)){
          GPLID =  Ubc_GPLs[u]
          df = get(GPLID) 
          which(rownames(df) == GEN)
          df_aux = df[GEN,]
          # Convertimos el df en un vector
          v_aux = as.numeric(as.vector(df[GEN,]))
          GEN_stat_vector = c(GEN_stat_vector, v_aux)}
        # Extraidos los valores para un gen:
        # anotamos el nº de studios en los que aparece y calculamos el valor mediano
        # Nº de estudios en los que aparece = longitud vector almacenamiento menos los NAs
        nGSE_GEN = length(GEN_stat_vector) - sum(is.na(GEN_stat_vector))
        nGSEs = c(nGSEs,nGSE_GEN)
        # Obtenemos el valor mediano (sin contar NAs) y lo almacenamos:
        MEDIAN = median(GEN_stat_vector, na.rm = TRUE)
        MEDIAN_VALUE = c(MEDIAN_VALUE, MEDIAN)
      }
    }
    
    rm(df, df_aux, GEN, MEDIAN,nGSE_GEN, v_aux)

    # Construimos el df de salida:
    # SYMBOL_ID | MEDIAN_VALUE | nGSE
    naiveDF = data.frame(SYMBOL_ID,MEDIAN_VALUE,nGSEs)
    # Ordenamos por la mediana para generar el ranking:
    naiveDF = naiveDF[order(naiveDF$MEDIAN_VALUE),]
    naiveDF$RANKING <- seq.int(nrow(naiveDF))
    # Reordenamos las columnas
    naiveDF = naiveDF[,c(1,4,2,3)]
    
    #dir.create("naiveRank")
    #write.csv(naiveDF,paste0("./RESULTS_",stat,"/naiveRank_",stat,".csv")) 
    #write.csv(naiveDF,paste0("./naiveRank/naiveRank_",stat,".csv"))
    write.csv(naiveDF,paste0("./naiveRank/",SEX,"_naiveRank_",stat,".csv"))
  }
}  




## ¿DONDE CAEN NUESTROS GENES HKG?
# Para cada estadistico, leemos los naiveranking y extraemos las filas de los HKG
## EXPLORAMOS LOS RESULTADOS:
stats = c("CV","IQRmedian","MADmedian")
SEXOS = c("ALL", "MALE", "FEMALE")
stats = c("CV","IQRmedian","MADmedian")
for (i in 1:length(SEXOS)){
  SEX = SEXOS[i]
  for (j in 1:length(stats)){
    stat = stats[j]
    if (SEX == ""){
      rank <- read.csv(paste0("./naiveRank/naiveRank_",stat,".csv"), header = TRUE )
    }
    if (SEX != ""){
      rank <- read.csv(paste0("./naiveRank/",SEX,"_naiveRank_",stat,".csv"), header = TRUE )
    }
    text = paste0("\ncondicion: ",SEX," - ",stat,"\n______________________________________________________\n")
    cat(text)
    df_aux = rank %>% select(SYMBOL_ID, RANKING ,MEDIAN_VALUE, nGSEs) %>% filter(SYMBOL_ID %in% c("Hprt", "Gapdh","Ppia", "Ubc", "Rpl19", "Rn18s"))
    print(df_aux)
    write.csv(df_aux,paste0("./naiveRank/",SEX,"HKGnaiveRank",stat,".csv"),row.names = FALSE)
  }
}

# Las mini tablitas de salida las presentamos en slides del TFM resultados
#SEXHKGnaiveRankstat.csv

