#### LIBRERIAS
library (GEOquery); packageDescription ("GEOquery", fields = "Version") #"2.50.5"
library(Biobase)
###
library(affy) # Affymetrix
### Manejo data.frames 
library(dplyr)
library(readr)
library(sqldf)
library(sqldf)

### PROCESAMIENTO DIFERENCIAS EN SEXO ######
# Obj1. Identificar las GSMs correspondientes de cada sexo
# Obj2. Contar las muestras de cada sexo 
library (GEOmetadb); packageDescription ("GEOmetadb", fields = "Version") #"1.44.0"
###Get the database file: "GEOmetadb.sqlite"
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
getSQLiteFile()
file.info('GEOmetadb.sqlite')
dir ()
###Connect to the database - se crea la conexión (objeto conexión)
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')


## 1. LECTURA FICHERO PROCESSED GSEs - EXTRAEMOS AQUELLOS ESTUDIOS CON INFORMACION DE SEXO (26 GSEs)
# 1.1 Lectura fichero ProcessedGSEs.csv
library(readxl)
GSE_GPL_sexsamples <- read_excel("~/Documentos/TFM/CLUSTER/ProcessedGSEs/GSE_GPL_sexsamples.xlsx")

### TABLA nGSEs por GPL:
GPLID = "GPL16570"
query = paste0('SELECT COUNT(GSEID) FROM GSE_GPL_sexsamples WHERE InformacionSexo = "SI" AND GPL ="',GPLID,'"')
count_aux = sqldf(query) #df
print(count_aux)



# 1.2 Extraemos los GSEIDs de los estudios que incluyen la información del sexo:
query = paste0('SELECT GSEID, GPL, nSamples_adi, indicadorMT, muestrasSexo, anotacionSexo, sex_tags FROM GSE_GPL_sexsamples WHERE InformacionSexo = "SI"')
df_aux = sqldf(query) #df
GSEIDs = df_aux$GSEID  # character
GPLIDs = df_aux$GPL
indicadorMT = df_aux$indicadorMT
muestrasSexo = df_aux$muestrasSexo
anotacionSexo = df_aux$anotacionSexo
sex_tags = df_aux$sex_tags


## 2. IDENTIFICAR Y CONTAR EL SEXO DE LAS MUESTRAS
# 2.1 - Para cada estudio obtenemos los identificadores y la descripción y el sexo de sus muestras
### Creamos la tabla GSE | GSM | Sample_Source | characteristics
### En la columna de characteristics está la información del sexo de la muestra


GSMformat <- function(GSM_aux) {
  # Input salida de la query:
  # Output: lista GSM formateada para almacenar
  # Convertimos la lista en una cadena:
  q = GSM_aux
  GSM_string = ''
  for (i in 1:length(q)) {
    gsm = q[i]
    if (i == 1) { # sabemos que GPL570 no es la primera
      gsm = paste("\'",gsm,"\'", sep ='')  # necesitamos añadirle las comillas
      GSM_string = paste(GSM_string, gsm, sep = '' )
      }else {
        gsm = paste("\'",gsm,"\'", sep ='')  # necesitamos añadirle las comillas
        GSM_string = paste(GSM_string, gsm, sep = ', ' )
        }
  }
  return(GSM_string)
}

# Listas para almacenar la nueva informacion: femaleGSMs, maleGSMs, nSamples_Female, nSamples_Male
# que incorporaremos al df_aux del punto 1.2

femaleGSMs = c()
maleGSMs = c()
nFemales = c()
nMales = c()

## Para Hsa - sampleSex = "AMBOS/HOMBRES/MUJERES"
## Para Mmu - sampleSex = "AMBOS/MACHOS/HEMBRAS"
for (i in 1:length(GSEIDs)) {# GPL_Express/brain_GPLs
  print(i)
  GSE = GSEIDs[i]
  anotSexo = anotacionSexo[i]
  samplesSex = muestrasSexo[i]
  GPL = GPLIDs[i]
  
  if (anotSexo == "SI"){ # El sexo de la muestra está anotado en las caracteristicas de la muestra:
  sextag = sex_tags[i]
  }
  
  sample_query <- paste("SELECT DISTINCT gse.gse, gsm.gsm, gsm.gpl, gsm.source_name_ch1 as sample_source, gsm.characteristics_ch1 as characteristics ,gsm.organism_ch1 as organism",
                        " FROM",
                        " gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
                        " JOIN gse ON gse_gsm.gse=gse.gse",
                        " WHERE",
                        " gsm.source_name_ch1 LIKE '%adipo%' AND",
                        " gsm.gpl = '",GPL,"' AND",
                        " gse.gse ='",GSE,"'", sep ='')
  
  sampleinformation_df <- dbGetQuery(con, sample_query)
  
  if (samplesSex == "AMBOS") { 
    # El estudio incluye muestras de ambos sexos:
    # Extraemos las etiquetas del sexo de la muestra: Propia de cada estudio
    sextag_vector = unlist(strsplit(sextag,"///")) #  Separamos por el //
    male_sextag = sextag_vector[1]
    female_sextag = sextag_vector[2]
    
    get_male_samples <- paste('SELECT gsm, characteristics FROM sampleinformation_df WHERE',
                              ' sample_source LIKE "%adipo%" AND',
                              ' characteristics LIKE "%',male_sextag,'%" AND', # Incluimos todas las variantes de male sex_tags (dato conocido)
                              ' organism LIKE "%Mus musculus%"', sep='')
    
    male_samples = sqldf(get_male_samples)
    GSMaux = male_samples$gsm
    maleGSMs_string = GSMformat(GSMaux)
    maleGSMs = c(maleGSMs, maleGSMs_string)
    count_males = length(GSMaux)
    nMales = c(nMales, count_males)
    
    get_female_samples <- paste('SELECT gsm, characteristics FROM sampleinformation_df WHERE',
                                ' sample_source LIKE "%adipo%" AND',
                                ' characteristics LIKE "%',female_sextag,'%" AND', # Incluimos todas las variantes de female sex_tags (dato conocido)
                                ' organism LIKE "%Mus musculus%"', sep='')
    
    female_samples = sqldf(get_female_samples)
    GSMaux = female_samples$gsm # Extraemos las GSMs
    femaleGSMs_string = GSMformat(GSMaux)
    femaleGSMs = c(femaleGSMs, femaleGSMs_string) # Almacenamos las GSMs
    count_females = length(GSMaux) # Contamos las GSMs
    nFemales = c(nFemales, count_females) # Almacenamos el conteo
  }

  if (samplesSex == "MACHOS") {
    # Todas las muestras del estudio son HOMBRES - no es necesario el sex_tag
    GSMaux = sampleinformation_df$gsm
    maleGSMs_string = GSMformat(GSMaux)
    maleGSMs = c(maleGSMs, maleGSMs_string)
    count_males = length(GSMaux)
    nMales = c(nMales, count_males)
    
    # Estudio solo de hombres no tiene femaleGSMs:
    femaleGSMs = c(femaleGSMs, "NULL")
    nFemales = c(nFemales, 0)
  }
  
  if (samplesSex == "HEMBRAS") {
    # Todas las muestras del estudio son MUJERES - no es necesario el sex_tag
    GSMaux = sampleinformation_df$gsm
    femaleGSMs_string = GSMformat(GSMaux)
    femaleGSMs = c(femaleGSMs, femaleGSMs_string)
    count_females = length(GSMaux)
    nFemales = c(nFemales, count_females)
    
    # Estudio solo de mujeres no tiene maleGSMs:
    maleGSMs = c(maleGSMs, "NULL")
    nMales = c(nMales, 0)
  }
}

## RESULTADO: Df_aux
#[1] "GSEID"         "GPL"           "nSamples_adi"  "indicadorMT"   "muestrasSexo"  "anotacionSexo" "sex_tags"     

# Reordenamos las columnas
df_aux <- df_aux[c("GSEID","GPL","indicadorMT", "muestrasSexo", "anotacionSexo","sex_tags","nSamples_adi")]
# Añadimos la nueva informacion:
df_aux$nMales = nMales
df_aux$nFemales = nFemales

df_aux$MaleGSMs = maleGSMs
df_aux$FemaleGSMs = femaleGSMs

# Escribimos el resultado:
write.csv2(df_aux,file = 'GSE_Sexsamples.csv', row.names = FALSE)
























