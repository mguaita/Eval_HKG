#### LIBRERIAS
library (GEOquery); packageDescription ("GEOquery", fields = "Version") #"2.50.5"
library(Biobase)
###
library(affy) # Affymetrix
### Manejo data.frames 
library(dplyr)
library(readr)
library(sqldf)

################################################################################
#1.- INTRODUCCION DE VARIABLES AL ENTORNO:

# Plataforma que vamos a trabajar:
#GPLID ="GPL10558"

# Lectura del csv (tab separated) summary del paso 5
# Se procesa para obtener los GSEIDs de los datasets que se quieren descargar:
# estudios con al menos 10 muestras
# Separar los identificadores de los estudios monotejido de los multitejido en diferentes listas
# para procesarlos más tarde

### HUMANOS:
# Hsa - 
setwd(paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated/",GPLID))
if (GPLID == "GPL570"){
  summary <- read_delim("GPL570_summary.csv", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
}

if (GPLID == "GPL6244"){
  summary <- read_delim("GPL6244_summary.csv",
                        "\t", escape_double = FALSE, trim_ws = TRUE)
}

if (GPLID == "GPL10558" || GPLID == "GPL6947"){
  summary <- read_csv(paste0(GPLID,"_summary.csv"))
}

### RATONES:
# Mmu - setwd("/clinicfs/userhomes/mguaita/Mmu_adipose")
setwd("/clinicfs/userhomes/mguaita/Mmu_adipose")

GPLID = "GPL1261"
summary <- read_csv(paste0("./",GPLID,"/",GPLID,"_summary.csv"))

# Extraemos todos los identificadores de los estudios con al menos 10 muestras de tejido adiposo
query_allGSEIDs = paste('SELECT GSEIDs, nSamples_adi ,indicadorMT, adiGSMs FROM summary WHERE nSamples_adi >= 10', sep='')
df_aux = sqldf(query_allGSEIDs) #df
GSEIDs = df_aux$GSEIDs
indicadorMT = df_aux$indicadorMT
adiGSMs = df_aux$adiGSMs
write.csv(ProcessedGSEs, paste0("./",GPLID,"/ProcessedGSEs_",GPLID,".csv"), row.names = FALSE)

rm(df_aux, query_allGSEIDs, summary)

################################################################################
#2.-DESCARGA DE DATOS: FUNCION DE DESCARGA
getGEOData <- function(GSEID,GPLID){
  # GSEID - Identificador GSE de GEO, string
  # GPLID - Identificador de la plataforma de GEO, string
  
  gseName = GSEID
  
  # Descarga de datos
  # if getGPL = FALSE - no tenemos feature data
  gset <- getGEO(GSEID, GSEMatrix = TRUE)
  if (length(gset) > 1) idx <- grep(GPLID, attr(gset, "names")) else idx <- 1
  Eset <- gset[[idx]]
  
  assign(gseName, Eset)
  filename = paste(GSEID, ".Rdata", sep ="")
  save (list = gseName, file = paste("./",GPLID,"/GEORdata_",GPLID,"/", filename, sep = ""))
}

# Descarga de datos:
#dir.create(paste0("GEORdata_",GPLID))
GPLs = c("GPL1261", "GPL6246","GPL6887", "GPL6885", "GPL16570")
for (g in 1:length(GPLs)){
  GPLID = GPLs[g]
  ## Leemos el fichero de estudios a procesar:
  ProcessedGSEs <- read_csv(paste0("./",GPLID,"/ProcessedGSEs_",GPLID,".csv"))
  GSEIDs = ProcessedGSEs$GSEIDs
  for (i in 1:length(GSEIDs)) { # puede haber un fallo de conexión
    cat(i, " : ", GSEIDs[i], "\n")
    tryCatch(getGEOData(GSEIDs[i],GPLID), silent = FALSE, 
             error = function(e) {
               n_Errores = n_Errores + 1
               append(index_error, GSEIDs[i-1])
               cat("Error: ", GSEIDs[i-1])}) # Mostramos el ID anterior
  }
}

################################################################################


################################################################################
#3.- CARGA DE LOS DATOS: Rdata
GPLID ='GPL570'
# Hsa:
path_ProcessedGSEs = paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated/",GPLID,"/ProcessedGSEs_",GPLID,".csv")
# Mmu:
path_ProcessedGSEs = paste0("/clinicfs/userhomes/mguaita/Mmu_adipose/",GPLID,"/ProcessedGSEs_",GPLID,".csv")

ProcessedGSEs <- read_csv(path_ProcessedGSEs)
GSEIDs=ProcessedGSEs$GSEIDs
indicadorMT = ProcessedGSEs$indicadorMT
adiGSMs = ProcessedGSEs$adiGSMs

## IMPORTANTE: 
## revisar que le primer estudio está completo - si no es así, modificar el 
## el orden de las filas en PorcessedGSEs
## Revisar que los indices de las listas coinciden con el contenido del fichero

# Hsa:
setwd(paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated/",GPLID,"/GEORdata_",GPLID))
# Mmu:
#setwd(paste0("/clinicfs/userhomes/mguaita/Mmu_adipose/",GPLID,"/GEORdata_",GPLID))
wd = getwd()
for (file in list.files(wd)) {load(file)}

################################################################################


################################################################################
# 4.- PROCESAMIENTO DE EXPRESSIONSET-MATRIZ DE EXPRESION --> C.V. DATAFRAME
################################################################################
### FUNCIONES DE ANOTACION ESPECIFICAS DE CADA PLATAFORMA:
## Si se quiere incorporar más plataformas al estudio, solo es necesario 
## incorporar su función de anotación especifica y conocer el nº máximo de genes
## que podemos anotar/identificar.
################################################################################
# Toma como argumentos:
# - objeto expressionSet, ya que extraemos los identificadores SYMBOL del fdata
# - GSEMat, matriz con los valores de expresión de las sondas que vamos a anotar
#
# Procesamiento propio de cada plataforma para procesar la anotación + control de nulos
# Incorporamos los identificadores SYMBOLs a la matriz de expresión
#
# Devolvemos la matriz anotada GSEMat
## ANOTACION HUMANOS:
anotacion_gpl570 <- function(Eset, GSEMat){
  ### CORREGIR los SYMBOLs con sinonimos
  # Mostramos la anotación disponible:	
  # print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Gene Symbol" # puede contener sinonimos
  
  ## ANOTACION SONDAS CONTROL:
  # PROBEIDs = rownames(GSEMat)
  # RNA18S_control_probes = c("AFFX-HUMRGE/M10098_3_at","AFFX-HUMRGE/M10098_5_at","AFFX-HUMRGE/M10098_M_at")
  # for (c in 1:length(RNA18S_control_probes)){
  #   cProbe = RNA18S_control_probes[c]
  #   # Anotamos el IDENTIFICADOR DEL RNA18S
  #   SYMBOLs[which(PROBEIDs == cProbe)] = "RNA18S"
  # }
  # 
  # Añadimos los identificadores a la matriz de expresión genica
  #GSEMat$SYMBOL = SYMBOLs # puede contener sinonimos
  # Antes de añadir la columna a la matriz nos quedamos con un identificador unico:
  processedSYMBOLs = c() # lista para almacenar los identificadores unicos
  for (j in 1:length(SYMBOLs)) {
    ID = SYMBOLs[j]
    if (grepl("///",ID)){
      ID = gsub(' ','',ID) # eliminamos espacios
      ID_vector = unlist(strsplit(ID,"///")) #  Separamos por el //
      Symbol = ID_vector[1]
    }else{
      Symbol = SYMBOLs[j]
    }
    processedSYMBOLs = append(processedSYMBOLs, Symbol)
  }
  # Añadimos los identificadores a la matriz:
  GSEMat$SYMBOLS = processedSYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 8890
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45785    15 
  
  return(GSEMat)
}

anotacion_gpl6244 <- function(Eset, GSEMat){
  # Mostramos la anotación disponible:	
  #print(colnames(fData(Eset)))
  #SYMBOLs = fData(Eset)$"SYMBOL"
  #SYMBOLs_bck = SYMBOLs
  
  gene_assignments = fData(Eset)$"gene_assignment"
  SYMBOLs = c()
  for (j in 1:length(gene_assignments)) {
    GA = gene_assignments[j]
    if (GA != "---"){
      GA = gsub(' ','',GA) # eliminamos espacios
      GA_vector = unlist(strsplit(GA,"//")) #  Separamos por el //
      Symbol = GA_vector[2]
    }else{
      Symbol = "NA"
    }
    SYMBOLs = append(SYMBOLs, Symbol)
  }
  
  GSEMat$SYMBOL = SYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "NA")) # 8197 
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != "NA"
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
  
}

anotacion_gpl10558 <- function(Eset, GSEMat){# Mostramos la anotación disponible:	
  # print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Symbol"
  
  # # Corregimos anotación rRNA 18S - LOC100008588 (ILMN_3243593)
  RNA18S_LOC = 'LOC100008588'
  if (is.element(RNA18S_LOC, SYMBOLs) == TRUE){ # Si esta la sonda, anotamos
     SYMBOLs[match(RNA18S_LOC, SYMBOLs)] = 'RNA18S5'} # Empleamos match porque solo hay una ocurrencia de la sonda, 1st match
   
  # Corregimos anotación de los meses:
  months = c("1-Dec","1-Mar","10-Mar","11-Mar","2-Mar","3-Mar","4-Mar","5-Mar","6-Mar","7-Mar","8-Mar","9-Mar")
  gen_correction = c("DEC1","MARCHF1","MARCHF10","MARCHF11","MARCHF2","MARCHF3","MARCHF4","MARCHF5","MARCHF6","MARCHF7","MARCHF8","MARCHF9")
  for (i in 1:length(months)) {
    month = months[i]
    gen_name = gen_correction[i]
    SYMBOLs[which(SYMBOLs == month)] = gen_name
  }
  
  ## Añadimos la anotación a la matriz:
  GSEMat$SYMBOL = SYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 33927
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
}

anotacion_gpl6947 <- function(Eset, GSEMat){
  # Mostramos la anotación disponible:	
  #print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Symbol"
  #SYMBOLs_bck = SYMBOLs
  
  ## Añadimos la anotación a la matriz:
  GSEMat$SYMBOL = SYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 33927
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
}

## ANOTACION RATONES:
anotacion_gpl1261 <- function(Eset, GSEMat){
  ### CORREGIR los SYMBOLs con sinonimos
  # Mostramos la anotación disponible:	
  # print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Gene Symbol" # puede contener sinonimos
  
  
  ## ANOTACION SONDAS CONTROL:
  PROBEIDs = rownames(GSEMat)
  # RNA18S_control_probes = c("AFFX-18SRNAMur/X00686_3_at","AFFX-18SRNAMur/X00686_5_at","AFFX-18SRNAMur/X00686_M_at")
  # for (c in 1:length(RNA18S_control_probes)){
  #   cProbe = RNA18S_control_probes[c]
  #   # Anotamos el IDENTIFICADOR DEL RNA18S
  #   SYMBOLs[which(PROBEIDs == cProbe)] = "Rn18s"
  # }
  
  ## ANOTACION Ppia:
  SYMBOLs[which(PROBEIDs == "1417451_a_at")] = "Ppia"
  # Añadimos los identificadores a la matriz de expresión genica
  #GSEMat$SYMBOL = SYMBOLs # puede contener sinonimos
  # Antes de añadir la columna a la matriz nos quedamos con un identificador unico:
  processedSYMBOLs = c() # lista para almacenar los identificadores unicos
  for (j in 1:length(SYMBOLs)) {
    ID = SYMBOLs[j]
    if (grepl("///",ID)){
      ID = gsub(' ','',ID) # eliminamos espacios
      ID_vector = unlist(strsplit(ID,"///")) #  Separamos por el //
      Symbol = ID_vector[1]
    }else{
      Symbol = SYMBOLs[j]
    }
    processedSYMBOLs = append(processedSYMBOLs, Symbol)
  }
  # Añadimos los identificadores a la matriz:
  GSEMat$SYMBOLS = processedSYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 5435
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
}

anotacion_gpl6246 <- function(Eset, GSEMat){
  # Mostramos la anotación disponible:	
  #print(colnames(fData(Eset)))
  #SYMBOLs = fData(Eset)$"SYMBOL"
  #SYMBOLs_bck = SYMBOLs
  
  gene_assignments = fData(Eset)$"gene_assignment"
  SYMBOLs = c()
  for (j in 1:length(gene_assignments)) {
    GA = gene_assignments[j]
    if (GA != "---"){
      GA = gsub(' ','',GA) # eliminamos espacios
      GA_vector = unlist(strsplit(GA,"//")) #  Separamos por el //
      Symbol = GA_vector[2]
    }else{
      Symbol = "NULL"
    }
    SYMBOLs = append(SYMBOLs, Symbol)
  }
  
  GSEMat$SYMBOL = SYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "NULL")) # 8197 
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != "NULL"
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 27359    49  
  
  return(GSEMat)
  
}

anotacion_gpl6887 <- function(Eset, GSEMat){
  ### CORREGIR los SYMBOLs con sinonimos
  # Mostramos la anotación disponible:	
  # print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Symbol" # puede contener sinonimos
  
  processedSYMBOLs = c() # lista para almacenar los identificadores unicos
  for (j in 1:length(SYMBOLs)) {
    ID = SYMBOLs[j]
    if (grepl(",",ID)){
      #print('SYMBOL corregido')
      ID = gsub(' ','',ID) # eliminamos espacios
      ID_vector = unlist(strsplit(ID,",")) #  Separamos por el //
      Symbol = ID_vector[1]
    }else{
      Symbol = SYMBOLs[j]
    }
    processedSYMBOLs = append(processedSYMBOLs, Symbol)
  }
  
  processedSYMBOLs[which(processedSYMBOLs == 'Hprt1')] = 'Hprt'
  
  # Añadimos los identificadores a la matriz de expresión genica
  GSEMat$SYMBOLS = processedSYMBOLs
  
  ## NO TIENE NULOS - Raro, los nulos los rellena con el codigo de ORF de la plataforma
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 5435
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
}

anotacion_gpl6885 <- function(Eset, GSEMat){
  ### CORREGIR los SYMBOLs con sinonimos
  # Mostramos la anotación disponible:	
  print(colnames(fData(Eset)))
  SYMBOLs = fData(Eset)$"Symbol" # puede contener sinonimos
  
  SYMBOLs[which(SYMBOLs == 'Gapd')] = 'Gapdh'
  # No tiene sinonimos lo9s alamcena en otra columna
  # Añadimos los identificadores a la matriz de expresión genica
  GSEMat$SYMBOLS = SYMBOLs
  
  ## NO TIENE NULOS - Raro, los nulos los rellena con el codigo de ORF de la plataforma
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "")) # 5435
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != ""
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 45782 13  
  
  return(GSEMat)
}

anotacion_gpl16570 <- function(Eset, GSEMat){
  # Mostramos la anotación disponible:	
  print(colnames(fData(Eset)))
  #SYMBOLs = fData(Eset)$"SYMBOL"
  #SYMBOLs_bck = SYMBOLs
  
  gene_assignments = fData(Eset)$"gene_assignment"
  SYMBOLs = c()
  for (j in 1:length(gene_assignments)) {
    GA = gene_assignments[j]
    if (GA != "---"){
      GA = gsub(' ','',GA) # eliminamos espacios
      GA_vector = unlist(strsplit(GA,"//")) #  Separamos por el //
      Symbol = GA_vector[2]
    }else{
      Symbol = "NULL"
    }
    SYMBOLs = append(SYMBOLs, Symbol)
  }
  
  GSEMat$SYMBOL = SYMBOLs
  
  # Sabemos de la existencia de NAs
  table((GSEMat$SYMBOL == "NULL")) # 12100 
  
  # Los eliminamos:
  filt_NA = GSEMat$SYMBOL != "NULL"
  GSEMat = GSEMat[filt_NA,]
  dim(GSEMat) # 27359    49  
  
  return(GSEMat)
  
}
################################################################################
### FUNCIONES DEL FLUJO DE PROCESAMIENTO DE DATOS:
################################################################################
# ANOTACION GEO (fData)
get_geneExpression <- function(GPLID, GSEID, Eset, MT, adipose_samples = "NA") {
  # Función que nos devuelve los valores de expresión medianos para un gen
  # Sustituimos el PROBEID de las sondas por el gen al que se corresponden
  # Colapasamos GSEMat agrupando por GENE SYMBOL, nos quedamos con la mediana como valor de expresión
  
  GSEMat = as.data.frame(exprs(Eset)) # 47323
  
  ## Procesamos la GSEMat para quedarnos con las muestras correspondientes al tejido adiposo:
  # Si el indicadorMT es multitejido procesamos la matriz de expresión genica para quedarnos
  # con las columnas de tejido adiposo:
  if (MT == 'MULTITEJIDO'){ # Extraemos las columnas de las muestras de tej.adiposo
    print("Getting adipose samples from GSEMat")
    # Con los identificadores de las muestras de tejido adiposo, estan en formato string
    # Lo covertimos a vector para seleccionar las columnas del df
    # Extraemos las columnas de tejido adiposo
    GSEMat = GSEMat[, adipose_samples]
  }
  
  ### CORREGIR MATRICES CON VALORES NEGATIVOS
  GSEMat = negativevalues_correction(GSEMat)
  
  PROBEIDs = rownames(GSEMat)
  
  ### FUNCIONES DE ANOTACION: Extraemos los SYMBOL del fdata(Eset) y los incorporamos a la matriz de expresion:
  
  ## HUMANOS:
  if (GPLID == "GPL570"){GSEMat = anotacion_gpl570(Eset, GSEMat)}
  if (GPLID == "GPL6244"){GSEMat = anotacion_gpl6244(Eset, GSEMat)}
  if (GPLID == "GPL10558"){GSEMat = anotacion_gpl10558(Eset, GSEMat)}
  if (GPLID == "GPL6947"){GSEMat = anotacion_gpl6947(Eset, GSEMat)}
  
  ## RATONES:
  if (GPLID == "GPL1261"){GSEMat = anotacion_gpl1261(Eset, GSEMat)}
  if (GPLID == "GPL6246"){GSEMat = anotacion_gpl6246(Eset, GSEMat)}
  if (GPLID == "GPL6887"){GSEMat = anotacion_gpl6887(Eset, GSEMat)}
  if (GPLID == "GPL6885"){GSEMat = anotacion_gpl6885(Eset, GSEMat)}
  if (GPLID == "GPL16570"){GSEMat = anotacion_gpl16570(Eset, GSEMat)}
  
  print("Annotated GSEMat")
  ###
  
  # Procesamos la GSEMat -> data_clean (es un df con los genes y sus valores de expresión por muestra
  n = unique(GSEMat$SYMBOL)
  cat("Nº identificadores SYMBOL:",length(n),"\n") # 21495
  n_samples = (dim(GSEMat))[2] - 1 # -1 porque la ultima columna es el identificador symbol
  cat("Nº de muestras:", n_samples, "\n")
  
  # Colapasamos GSEMat agrupando por GENE SYMBOL, nos quedamos con la mediana:
  print("Aggregating expression values by SYMBOL")
  geneExpressionDF = aggregatebySYMBOL(GSEMat, n_samples)
  
  ## Guardamos la matriz de expresión génica anotada:
  print("Writing GSEdf")
  writeGeneExpressiondata(GPLID,geneExpressionDF,GSEID)
  
  return(geneExpressionDF)
}

add_resPAR_to_finalRes <- function(final_res, res_PAR) {
  # Función para añadir la columna del estudio procesado al df de resultados
  ## PODEMOS EMPLEAR LA FUNCION cbind porque comparten el orden de los identificadores SYMBOL
  # table(GSE79434[,1] == GSE37514[,1])
  # Unimos los df por columnas con cbind
  final_res = cbind(final_res, res_PAR)
  return (final_res)
}

getwrite_cleanGeneExpressionDF <- function(GPLID, GSEIDs, indicadorMT, adiGSMs){
  # Procesa el GSE R data de cada estudio , anota los genes, calcula la expresión génica mediana de cada gen 
  # y escribe el dataframe de cada estudio.
  # Llama a las funciones 1.get_geneExpression y 2.writeGeneExpressiondata (anidada en la 1ra)
  
  ## Obtenemos las matrices de expresión génica ya anotadas
  nMAXGENES = 0
  # COMENZAMOS BUCLE:
  for (i in 1:length(GSEIDs)) {
    GSEID = GSEIDs[i]
    cat(i,":", GSEID, "\n")
    # Procesamiento diferencial estudios monotejido/multitejido::
    ## Procesamos la información del estudio para quedarnos con las muestras correspondientes al tejido adiposo:
    # Si el indicadorMT es multitejido procesamos la matriz de expresión genica para quedarnos
    # con las columnas de tejido adiposo:
    if (indicadorMT[i] == 'MULTITEJIDO'){ 
      print("Estudio MULTITEJIDO")
      # Extraemos las columnas de las muestras de tej.adiposo
      # Extraemos los identificadores de las muestras de tejido adiposo, estan en formato string
      # Lo covertimos a vector para seleccionar las columnas del df
      GSM_IDs = adiGSMs[i]
      GSM_IDs = gsub(' ','',GSM_IDs)
      GSM_IDs = gsub("\'",'',GSM_IDs)
      adipose_samples = unlist(strsplit(GSM_IDs,","))
      
      # Obtención de la expresión génica mediana/Anotacion SYMBOL
      geneExpressDF = get_geneExpression(GPLID, GSEIDs[i], get(GSEIDs[i]), indicadorMT[i], adipose_samples)
    }else
    {# adipose_samples to NULL
      # Obtención de la expresión génica mediana/Anotacion SYMBOL
      geneExpressDF = get_geneExpression(GPLID, GSEIDs[i], get(GSEIDs[i]), indicadorMT[i])
    }
    ## Obtenemos el nº maximo de genes unicos que podemos identificar:
    nGenes = dim(geneExpressDF)[1]
    if(nGenes > nMAXGENES){
      GSE_nmax = GSEID
      nMAXGENES = nGenes
    }
  }
  cat("Nº de genes unicos identificados:",nMAXGENES,"\n",
      "Estudio con más genes identificados:", GSE_nmax,"\n",
      "Gen expression data anotation done. See GSEdf directory.")
  
}

geneExpressionDF_reader <- function(GPLID, GSEID){
  # Lectura de los dfs ya anotados
  # Necesitamos conocer la estructura del fichero y que se utiliza correctamente
  ## GPLID - para construir la ruta
  ## GSEID - para leer el fichero del estudio
  # Carpeta GEORdata_GPLID/GSEdf/GSEID.csv
  geneExpressionDF <- read.csv(paste0(GPLID,"/GEORdata_",GPLID,"/GSEdf_",GPLID,"/",GSEID,".csv"), header = TRUE, row.names = 1)
  return(geneExpressionDF)
} 

geneExpressionDF_loopprocessing <- function(GPLID, GSEIDs, PAR, nMAXGENES){
  # Procesamiento de matrices ya anotadas
  # COMENZAMOS BUCLE:
  for (i in 1:length(GSEIDs)) {
    GSEID = GSEIDs[i]
    cat(i,':',GSEID,'\n')
    ## Leemos el fichero: llama al reader - necesita saber la GPL
    geneExpressDF = geneExpressionDF_reader(GPLID, GSEID)
    
    # n_muestras - necesario para delimitar el número de cols sobre las que calcular los estadisticos
    n_samples = (dim(geneExpressDF))[2]
    
    # PROCESAMOS SEGUN EL PARAMETRO QUE QUERAMOS CALCULAR:
    if (PAR == "CV"){res_PAR = get_CV(GSEID, geneExpressDF, n_samples)}
    
    if (PAR == "IQRmedian"){res_PAR = get_IQRmedian(GSEID, geneExpressDF, n_samples)}
    
    if (PAR == "MADmedian"){res_PAR = get_MADmedian(GSEID, geneExpressDF, n_samples)}
    
    # Sabemos que tenemos: 54675 sondas y 23520 genes anotados (con la info de GEO)
    # Sabemos que tenemos 21194 genes anotados con el paquete de Bioconductor
    
    # PRIMERA ITERACION:
    if (i == 1) { 
      # Iniciamos el df final resultado: (hacemos copia solo si es correcto)
      final_res = res_PAR			
      if(dim(res_PAR)[1] == nMAXGENES) { # Data frame normativo
        res_bckp = res_PAR # para luego hacer "outer join" con los dfs incompletos -> genes sin datos devuelven un NA
      }
    }
    # Puede darse el caso que el 1er df no sea normativo:
    
    # SEGUNDA ITERACION
    ## ES MAS SENCILLO ORDENAR EL FICHERO ProcessedGSEs para que el 1er estudio procesado, sea
    ## un estudio completo, mantenemos el codigo por si acaso pero no deberia ejecutarse este fragmento:
    # A.- PRIMER ESTUDIO ES NORMAL: 
    #if (i == 2 && dim(final_res)[1] == nMAXGENES) { # JUNTAMOS TAL CUAL
    #  final_res = add_resPAR_to_finalRes(final_res, res_PAR)
    #}
    
    # B.- PRIMER ESTUDIO ES UN CASO ESPECIAL: 
    #if (i == 2 && dim(final_res)[1] != nMAXGENES) { # El primer df no tiene la longitud que toca, lo juntamos tal cual
    ## Guardamos el 2do df, que es normal, como back up
    # para luego hacer "outer join" con los dfs incompletos -> genes sin datos devuelven un NA
    # res_bckp = res_PAR 
    # Hacemos el join con un df normal para completar la información faltante con valores nulos NAs
    # Por la función el df resultante tiene los genes en los nombres de las columnas
    # con TRUE forzamos a incluir todas las filas
    #final_res <- merge(res_PAR, final_res, by= "row.names", all= TRUE) 
    
    # Cambiamos el orden: GSE2 (sin NAs) | GSE1 (con NAs)
    #}
    
    # SIGUIENTES ITERACIONES:
    if (i >= 2 && dim(res_PAR)[1] == nMAXGENES) { # Normativo, lo adjuntamos normal
      ## Podemos unirlos directamente porque tienen la misma dimensión y el mismo orden de genes
      final_res = add_resPAR_to_finalRes(final_res, res_PAR)
      print("Columna añadida al df final.")
    }
    
    if (i >= 2 && dim(res_PAR)[1] < nMAXGENES) { # Tiene menos genes - requiere más pre-procesamiento para unirlo a res
      # Hacemos el join con un df normal para completar la información faltante con valores nulos NAs	
      # Por la función el df resultante tiene los genes en los nombres de las columnas
      # con TRUE forzamos a incluir todas las filas
      join_df <- merge(res_bckp, res_PAR, by= "row.names", all= TRUE) 
      
      # obtenemos un df de 3 columnas: 1. Gene Symbol, 2. res_bckp$GSE.x , 3. res_PAR$GSE.y
      # solo nos interesa la columna del estudio excepcional: la tercera - GSEaccesion.y
      join_df_3 = join_df %>% select((3))
      colnames(join_df_3) = GSEID # quitamos el sufijo ".y"
      
      # Los juntamos:
      final_res = add_resPAR_to_finalRes(final_res, join_df_3)
      print("Columna añadida al df final.")
    }
    
    # Mostramos el estado actual del resultado:
    cat("Columnas df final:", colnames(final_res), "\n")
    
  }
  
  ### ELIMINAR LA COLUMNA ROW.NAMES
  colnames(final_res)
  final_res = final_res[ , !names(final_res) %in% c("Row.names")]
  
  ## BUCLE FINALIZADO:
  print(dim(final_res)) 
  
  # Escribimos el resultado final:
  writeOutput(final_res, GPLID, PAR)
  
  return(final_res)
  
}


################################################################################
### FUNCIONES DE PROCESAMIENTO MATEMATICO:
################################################################################
negativevalues_correction <- function(GSEMat){
  # Funcion para corregir las medidas de expresión negativas.
  if (min(GSEMat, na.rm = TRUE) < 0) { # Si la matriz tiene valores negativos
    minimum_value = min(GSEMat, na.rm = TRUE) # obtenemos el minimo
    GSEMat = GSEMat + abs(minimum_value) # se los sumamos a la matriz
  }
  return(GSEMat)
}

aggregatebySYMBOL <- function(GSEMat, n_samples){
  ## Funcion para agregar las medidas de las sondas que se corresponden con el mismo gen
  ## Nos quedamos con el valor mediano de la expresión medida.
  ## Empleamos el dato n_samples para excluir la última columna que son los SYMBOL
  # Colapasamos GSEMat agrupando por GENE SYMBOL, nos quedamos con la mediana:
  geneExpressionDF = aggregate(GSEMat[,1:n_samples], list(SYMBOL = GSEMat$SYMBOL), median)
  rownames(geneExpressionDF) = geneExpressionDF$SYMBOL
  
  geneExpressionDF = geneExpressionDF[,-1] # eliminamos la primera columna del identificador gene symbol
  print(dim(geneExpressionDF))
  return(geneExpressionDF)
}

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
  # La ultima columna se corresponde con el estadistico que queremos calcular:
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
  d$MADmedian = apply(d[,1:n_samples], 1, function(x){ 1.4826*(mad(x, na.rm = TRUE)/median(x, na.rm = TRUE))})
  
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

writeOutput <- function(final_res, GPLID, PAR) {
  # Escribimos en su propia carpeta:
  dir.create(paste0(GPLID,"/GEORdata_",GPLID,"/Result"))
  write.csv(final_res, paste0("./",GPLID,"/GEORdata_",GPLID,"/Result/",GPLID,"_",PAR,".csv"), row.names = TRUE)
}
################################################################################


####################################
########   MAIN PROGRAMME ##########
####################################
# Para procesar todos los datos de forma automatica estandarizamos el nombre de los directorios 
# EL directorio de trabajo es "Automated" contienen una carpeta por plataforma "GPLID"
# Y los Rdata de lso estudios estan en GPLID/GEORdata_GPLID

GPLID = "GPL570"
# Definimos el numero maximo de genes que somos capaces de identificar:
nMAXGENES = 22881
# Hsa:
setwd(paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated"))
# Mmu:
#setwd("/clinicfs/userhomes/mguaita/Mmu_adipose")
#1. OBTENEMOS LAS MATRICES DE EXPRESION GENICA ANOTADAS
## Procesamiento del objeto GSE Rdata, anota los genes, calcula la expresión génica mediana de cada gen
## Llama a las funciones 1.get_geneExpression y 2.writeGeneExpressiondata (anidada en la 1ra)
## Obtenemos las matrices de expresión génica ya anotadas
getwrite_cleanGeneExpressionDF(GPLID, GSEIDs, indicadorMT, adiGSMs)

#2. PROCESAMIENTO ESTADISTICO
# Definimos el parametro a procesar: "CV", "IQRmedian", "MADmedian"
# La función principal llama al reader de las matrices de expresión anotadas
parametros = c("CV", "IQRmedian", "MADmedian")
for (p in 1:length(parametros)) {
  PAR = parametros[p]
  print(PAR)
  final_res = geneExpressionDF_loopprocessing(GPLID, GSEIDs, PAR, nMAXGENES)
}

rm(list = setdiff(ls(), lsf.str()))

rm(final_res)
rm(list=ls(pattern="GSE"))
rm(adiGSMs,file, GPLID, indicadorMT, nMAXGENES, p, PAR, parametros)



## VERSION AUTOMATIZADA COMPLETA: para procesar todos los estudios de todas las plataformas:
GPLs = c("GPL1261","GPL6246","GPL6887","GPL6885","GPL16570")
nMG = c(21495,24213,30866,18120,24647)

GPLs  = c("GPL6244", "GPL10558", "GPL6947")
nMG = c(23307,31426,25159)

for (g in 1:length(GPLs)){
  GPLID = GPLs[g]
  nMAXGENES = nMG[g]
  print(GPLID)
  
  #1.- LECTURA DE VARIABLES
  #path_ProcessedGSEs = paste0("/clinicfs/userhomes/mguaita/Mmu_adipose/",GPLID,"/ProcessedGSEs_",GPLID,".csv")
  path_ProcessedGSEs = paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated/",GPLID,"/ProcessedGSEs_",GPLID,".csv")
  print("Lectura fichero Processed_GSEs")
  ProcessedGSEs <- read_csv(path_ProcessedGSEs)
  GSEIDs=ProcessedGSEs$GSEIDs
  indicadorMT = ProcessedGSEs$indicadorMT
  adiGSMs = ProcessedGSEs$adiGSMs
  
  ## IMPORTANTE: 
  ## revisar que le primer estudio está completo - si no es así, modificar el 
  ## el orden de las filas en PorcessedGSEs
  ## Revisar que los indices de las listas coinciden con el contenido del fichero
  
  #2.- LECTURA DE DATOS
  print("Cargando GSEs R data")
  #setwd(paste0("/clinicfs/userhomes/mguaita/Mmu_adipose/",GPLID,"/GEORdata_",GPLID))
  setwd(paste0("/clinicfs/userhomes/mguaita/Human/adipose/Automated/",GPLID,"/GEORdata_",GPLID))
  wd = getwd()
  for (file in list.files(wd)) {load(file)}
  
  
  #setwd("/clinicfs/userhomes/mguaita/Mmu_adipose")
  setwd("/clinicfs/userhomes/mguaita/Human/adipose/Automated")
  getwrite_cleanGeneExpressionDF(GPLID, GSEIDs, indicadorMT, adiGSMs)

  #3. PROCESAMIENTO ESTADISTICO
  # Definimos el parametro a procesar: "CV", "IQRmedian", "MADmedian"
  # La función principal llama al reader de las matrices de expresión anotadas
  parametros = c("CV", "IQRmedian", "MADmedian")
  for (p in 1:length(parametros)) {
    PAR = parametros[p]
    print(PAR)
    final_res = geneExpressionDF_loopprocessing(GPLID, GSEIDs, PAR, nMAXGENES)
  }
  ## Limpiamos el escritorio de objetos y variables, solo dejamos las funciones:
  rm(final_res)
  rm(list=ls(pattern="GSE"))
  rm(adiGSMs,file, GPLID, indicadorMT, nMAXGENES, p, PAR, parametros)
  
}

