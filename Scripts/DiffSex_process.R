###################################
#   MASTER SCRIPT ANALISIS SEXOS  #
###################################
#### LIBRERIAS
library (GEOquery); packageDescription ("GEOquery", fields = "Version") #"2.50.5"
library(Biobase)
###
library(affy) # Affymetrix
### Manejo data.frames 
library(dplyr)
library(readr)
library(sqldf)
library(stringr)

###########################################################
# Adaptación del MasterScript para realizar el
# análisis basado en sexos.
##########################################################


################################################################################
### FUNCIONES DEL FLUJO DE PROCESAMIENTO DE DATOS:
################################################################################

geneExpressionDF_reader <- function(GPLID, GSEID){
  # Lectura de los dfs ya anotados
  # Necesitamos conocer la estructura del fichero y que se utiliza correctamente
  ## GPLID - para construir la ruta
  ## GSEID - para leer el fichero del estudio
  # Carpeta GEORdata_GPLID/GSEdf/GSEID.csv
  geneExpressionDF <- read.csv(paste0(GPLID,"/GEORdata_",GPLID,"/GSEdf_",GPLID,"/",GSEID,".csv"), header = TRUE, row.names = 1)
  return(geneExpressionDF)
} 

get_SexedSamples <- function(GSEdf, GSMs){
  # toma la matriz de expresión génica del estudio y con los identificadores GSM
  # extrae las muestra del sexo de interés
  GSM_IDs = gsub(' ','',GSMs)
  GSM_IDs = gsub("\'",'',GSM_IDs)
  sex_samples = unlist(strsplit(GSM_IDs,","))
  
  sexedGSEdf = GSEdf[, sex_samples]
  
  return(sexedGSEdf)
}

get_RESPAR <- function(PAR, GSEID,GSEdf, n_samples){
  # Dirige el flujo de trabajo
  if (PAR == "CV"){res_PAR = get_CV(GSEID, GSEdf, n_samples)}
  
  if (PAR == "IQRmedian"){res_PAR = get_IQRmedian(GSEID, GSEdf, n_samples)}
  
  if (PAR == "MADmedian"){res_PAR = get_MADmedian(GSEID, GSEdf, n_samples)}
  
  return(res_PAR)
}

add_resPAR_to_finalRes <- function(final_res, res_PAR) {
  # Función para añadir la columna del estudio procesado al df de resultados
  
  # Unimos los df por columnas con cbind
  final_res = cbind(final_res, res_PAR)
  return (final_res)
}


################################################################################
### FUNCIONES DE PROCESAMIENTO MATEMATICO:
################################################################################

get_CV <- function(GSEID, d, n_samples) { # dos casos A) dim(GSEMat)[1] == n_sondas, B) dim(GSEMat)[1] != n_sondas
  # Toma el dataframe de expresión génica  - las filas son genes, columnas las muestras
  # calcula la sd, media y el C.V por filas
  
  ## PROCESAMIENTO PARA EL COEFICIENTE DE VARIACION: 
  # n_samples para no incluir en el procesamiento los valores añadidos
  ## IGNORAMOS LOS NAs
  d$sd = apply(d, 1, sd, na.rm = TRUE) # 1-filas
  d$mean = apply(d[,1:n_samples], 1, mean, na.rm = TRUE)
  d$C.V = apply(d[,1:n_samples], 1, function(x){sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)})
  
  #writeGSEExp_sd_mean(d,GSEID)
  
  ## EXTRAEMOS EL C.V.:
  # hacemos select de la última columna, n_samples+3 pq añadimos 3 columnas: sd-media-C.V.
  res_aux = d %>% select((n_samples+3)) 
  colnames(res_aux) = GSEID # 'restauramos' el nombre de la columna C.V por el GSEID del estudio 
  
  return(res_aux)
}

get_IQRmedian <- function(GSEID, d, n_samples){
  d$IQR = apply(d, 1, IQR, na.rm = TRUE) # 1-filas
  d$median = apply(d[,1:n_samples], 1, median, na.rm = TRUE)
  d$medianIQR = apply(d[,1:n_samples], 1, function(x){0.75*(IQR(x, na.rm = TRUE)/median(x, na.rm = TRUE))})
  
  ## EXTRAEMOS EL IQR/median:
  # hacemos select de la última columna, n_samples+3 pq añadimos 3 columnas: IQR-median-medianIQR.
  res_aux = d %>% select((n_samples+3)) 
  colnames(res_aux) = GSEID # 'restauramos' el nombre de la columna C.V por el GSEID del estudio 
  
  return(res_aux)
}

get_MADmedian <- function(GSEID, d, n_samples){
  d$MAD = apply(d, 1, mad, na.rm = TRUE) # 1-filas
  d$median = apply(d[,1:n_samples], 1, median, na.rm = TRUE)
  d$MADmedian = apply(d[,1:n_samples], 1, function(x){1.4826*(mad(x, na.rm = TRUE)/median(x, na.rm = TRUE))})
  
  ## EXTRAEMOS EL MAD/median:
  # hacemos select de la última columna, n_samples+3 pq añadimos 3 columnas: IQR-median-medianIQR.
  res_aux = d %>% select((n_samples+3)) 
  colnames(res_aux) = GSEID # 'restauramos' el nombre de la columna C.V por el GSEID del estudio 
  
  return(res_aux)
}

################################################################################
### FUNCIONES DE ESCRITURA DE DATOS
################################################################################
writeGSEExp_sd_mean <- function(df,GSEID) {
  dir.create("GSEdf_statistics")
  filename = paste0("./GSEdf_statistics/", GSEID, ".csv")
  write.csv(df, filename, row.names = TRUE)
}

writeGeneExpressiondata <- function(GPLID, df,GSEID) {
  dir.create(paste0(GPLID,"/GEORdata_",GPLID,"/GSEdf_",GPLID))
  filename = paste0("./",GPLID,"/GEORdata_",GPLID,"/GSEdf_",GPLID,"/",GSEID,".csv")
  write.csv(df, filename, row.names = TRUE)
}

writeOutput <- function(final_res, GPLID, PAR, SEX) {
  dir.create(paste0(GPLID,"/GEORdata_",GPLID,"/Result"))
  write.csv(final_res, paste0("./",GPLID,"/GEORdata_",GPLID,"/Result/",SEX,"_",GPLID,"_",PAR,".csv"), row.names = TRUE)
}
################################################################################



#1.- INTRODUCCION DE VARIABLES AL ENTORNO:

# Lectura del csv (SEMICOLON separated) 
library(readr)
GSE_Sexsamples <- read_csv("GSE_Sexsamples.csv")

# Extraemos todos los identificadores de los estudios con al menos 10 muestras de tejido adiposo
### TABLA nGSEs por GPL:
GPLID = "GPL570"
query = paste0('SELECT * FROM GSE_Sexsamples WHERE GPL ="',GPLID,'"')
#query = paste0('SELECT nSamples_adi, nMales, nFemales FROM GSE_Sexsamples WHERE GPL ="',GPLID,'"')
df_aux = sqldf(query) #df
#a =  lapply(count_aux, sum, na.rm = TRUE)

GSEIDs = df_aux$GSEID
muestrasSexo = df_aux$muestrasSexo
maleGSMs = df_aux$MaleGSMs
femaleGSMs = df_aux$FemaleGSMs

rm(df_aux, query)


####################################
########   MAIN PROGRAMME ##########
####################################
#Hsa -
GPLs = c("GPL570", "GPL6244","GPL10558","GPL6947")
NMGENES = c(22880,23307,31426,25159)
#Mmu -
GPLs = c("GPL1261","GPL6246","GPL6887","GPL6885","GPL16570")
NMGENES = c(21496,24213,30866,18120,24647)

SEXOS = c("MALE", "FEMALE")

for (s in 1:length(SEXOS)){ # Ejecutamos el bucle completo 2 veces, una por cada sexo
  print(SEXOS[s])
  for (g in 1:length(GPLs)) { # Ejecutamos este bucle 4 veces, una por plataform
    GPLID = GPLs[g]
    nMAXGENES = NMGENES[g]
    print(GPLID)
    
    query = paste0('SELECT * FROM GSE_Sexsamples WHERE GPL ="',GPLID,'"')
    df_aux = sqldf(query)
    GSEIDs = df_aux$GSEID
    cat("Nº de estudios:",length(GSEIDs),"\n")
    muestrasSexo = df_aux$muestrasSexo
    maleGSMs = df_aux$MaleGSMs
    femaleGSMs = df_aux$FemaleGSMs
    
    # BUCLE PARA PROCESAR LOS 3 PARAMETROS:
    parametros = c("CV", "IQRmedian", "MADmedian")
    
    for (p in 1:length(parametros)) { # Ejecutamos este bucle 3 veces, una por cada estadistico
      PAR = parametros[p]
      print(PAR)
      
      # Iniciamos la variable contador - cuenta los estudios que se han incluido en el analisis
      # Para dirigir el loop de construccion del DF final:
      n_ProcessedGSEs = 0 
      
      for (i in 1:length(GSEIDs)) { # Ejectutamos este bucle tantas veces como estudios se interrogran
        GSEID = GSEIDs[i]
        cat(i,".",GSEID,"\n")
        ## Seleccionamos las etiquetas identificadoras de las muestras:
        if (SEXOS[s] == "MALE"){sexedGSMs = maleGSMs[i]}
        if (SEXOS[s] == "FEMALE"){sexedGSMs = femaleGSMs[i]}
        
        ## Comenzamos el procesado de los estudios: pueden tener muestras del sexo que nos interesa o no,
        ## Si no las tienen entonces sexedGSMs = "NULL" - No procesamos el estudio
        ## Si sexedGSMs != NULL el estudio contiene muestras del sexo que estamos analizando:
        if (sexedGSMs != "NULL" && str_count(sexedGSMs, "GSM") > 1){
          ## Nos aseguramos de procesar estudios con más de una muestra
          # El estudio tiene muestras del sexo que estamos analizando, +1 al contador:
          n_ProcessedGSEs = n_ProcessedGSEs + 1
          
          # Leemos el df ya anotado del estudio:
          # Construir path del CLUSTER para la funcion READER
          GSEdf = geneExpressionDF_reader(GPLID, GSEID)
          
          # Extraemos las muestras
          sexedGSEdf = get_SexedSamples(GSEdf,sexedGSMs)
          
          n_samples = (dim(sexedGSEdf))[2]
          
          # GET RES_PAR -  PROCESAMOS SEGUN EL PARAMETRO QUE QUERAMOS CALCULAR:
          res_PAR = get_RESPAR(PAR, GSEID, sexedGSEdf, n_samples)
          
          # Almacenar FINAL_RES: CONDICIONALES DE CONSTRUCCION DEL RESULTADO FINAL
          ###################################################################################
          # PRIMERA ITERACION: Inicializamos el df final
          if (n_ProcessedGSEs == 1) { 
            # Iniciamos el df final resultado: (hacemos copia solo si es correcto)
            final_res = res_PAR
            if (SEXOS[s] == "MALE" && GPLID == "GPL10558"){
              ## CAso especial: ninguno de los 3/4 estudios con muestras de HOMBRE es normativo, revisión manual
              # el nº max (nMAXGENES) que identificamos en muestras de hombres en esta plataforma es: 31418
              # Cambiamos el valor de esta variable (= dim(GSEdf) del 2do estudio, 1ro de hombres)
              nMAXGENES = 31418
              # En mujerres no pasa porque el 1er estudio es normativo: nMAXGENES = 31426
            }
            
            if (dim(res_PAR)[1] == nMAXGENES) { # Data frame normativo
              res_bckp = res_PAR # para luego hacer "outer join" con los dfs incompletos -> genes sin datos devuelven un NA
            }
          }
          # Puede darse el caso que el 1er df no sea normativo:
          
          # SEGUNDA ITERACION
          # A.- PRIMER ESTUDIO ES NORMAL:
          if (n_ProcessedGSEs == 2 && dim(final_res)[1] == nMAXGENES && dim(res_PAR)[1]==dim(final_res)[1]) { # JUNTAMOS TAL CUAL
            final_res = add_resPAR_to_finalRes(final_res, res_PAR)
          }
          
          
          # B.- PRIMER ESTUDIO ES UN CASO ESPECIAL:
          if (n_ProcessedGSEs == 2 && dim(final_res)[1] != nMAXGENES) { # El primer df no tiene la longitud que toca, lo juntamos tal cual
            ## Guardamos el 2do df, que es normal, como back up
            # para luego hacer "outer join" con los dfs incompletos -> genes sin datos devuelven un NA
            res_bckp = res_PAR 
            # Hacemos el join con un df normal para completar la información faltante con valores nulos NAs
            # Por la función el df resultante tiene los genes en los nombres de las columnas
            # con TRUE forzamos a incluir todas las filas
            final_res <- merge(res_PAR, final_res, by= "row.names", all= TRUE) 
            final_res = final_res[ , !names(final_res) %in% c("Row.names")]
            
            # Cambiamos el orden: GSE2 (sin NAs) | GSE1 (con NAs)
          }
          
          
          # SIGUIENTES ITERACIONES:
          if (n_ProcessedGSEs > 2 && dim(res_PAR)[1] == nMAXGENES) { # Normativo, lo adjuntamos normal
            final_res = add_resPAR_to_finalRes(final_res, res_PAR)
          }
          
          if (n_ProcessedGSEs >= 2 && dim(res_PAR)[1] < nMAXGENES) { # Tiene menos genes - requiere más pre-procesamiento para unirlo a res
            # Hacemos el join con un df normal para completar la información faltante con valores nulos NAs	
            # Por la función el df resultante tiene los genes en los nombres de las columnas
            # con TRUE forzamos a incluir todas las filas
            join_df <- merge(res_bckp, res_PAR, by= "row.names", all= TRUE) 
            
            # Los juntamos:
            if (dim(final_res)[1] == dim(join_df)[1]){
              # obtenemos un df de 3 columnas: 1. Gene Symbol, 2. res_bckp$GSE.x , 3. res_PAR$GSE.y
              # solo nos interesa la columna del estudio excepcional: la tercera - GSEaccesion.y
              join_df_3 = join_df %>% select((3))
              colnames(join_df_3) = GSEID # quitamos el sufijo ".y"
              
              # Los juntamos:
              final_res = add_resPAR_to_finalRes(final_res, join_df_3)
              
            } 
            if (dim(final_res)[1] != dim(join_df)[1]){
              # Completamos con NAs tanto finale _res como RESPAR
              join_df <- merge(final_res, res_PAR, by= "row.names", all= TRUE)
              rownames(join_df) <- join_df$"Row.names"
              final_res = join_df
              final_res = final_res[ , !names(final_res) %in% c("Row.names")]
              if (dim(final_res)[1] > nMAXGENES){ nMAXGENES = dim(final_res)[1]}
            }
          }
          
          # Mostramos el estado actual del resultado:
          cat("res_colnames:", colnames(final_res), "\n")
          ###################################################################################
        }
      }
      
      # OBTENEMOS EL DF FINAL con los valores del estadistico de cada gen en cada estudio/sexo
      ### ELIMINAR LA COLUMNA ROW.NAMES
      colnames(final_res)
      # Si tiene columna "row.names" la eliminamos:
      #final_res = final_res[ , !names(final_res) %in% c("Row.names")]
      
      ## BUCLE FINALIZADO:
      print(dim(final_res)) 
      
      # Escribimos el resultado final:
      writeOutput(final_res, GPLID, PAR, SEXOS[s])
      
      ## Acabamos con un directorio "RESULT" con 3*4*2 = 24 ficheros, para cada sexo, las 4 GPLS y los 3 coeficientes
      # ./SEX_GPLID_PAR.csv
    }
  }
  
}





