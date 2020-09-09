######### PROCESAMIENTO POR TIPO DE PLATAFORMA ###########
# Obj. Identificar estudios monotejido y multitejido
# Obj2. Para los estudios multitejido, identificar las muestras correspondientes con tejido adiposo.
library(sqldf)

library (GEOmetadb); packageDescription ("GEOmetadb", fields = "Version") #"1.44.0"
###Get the database file: "GEOmetadb.sqlite"
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
getSQLiteFile()
file.info('GEOmetadb.sqlite')
dir ()
###Connect to the database - se crea la conexión (objeto conexión)
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')

############### PARTE 1 ###########################
# Procesamos por plataforma: 
## Seleccionamos de la base de datos relacional creada por el script GSE_GPL_relational.py
## gse_gpl_relational, muestra la relación estudios - plataforma
library(readr)
## Leemos la base de datos que hemos construido:
gse_gpl_relational <- read_delim("gse_gpl_relational.csv",
"\t", escape_double = FALSE, trim_ws = TRUE)

# Hsa - 
# Mmu - GPLs = c("GPL1261", "GPL6246", "GPL6887","GPL6885", "GPL16570")

# Creamos los directorios de trabajo:
for (i in 1:length(GPLs)){
  GPLID = GPLs[i]
  dir.create(GPLID)
}
## Interrogamos el TOP de plataformas de interés:
GPLID = 'GPL570'
query = paste('SELECT GSEID FROM gse_gpl_relational WHERE GPL = "',GPLID,'"', sep='')
df_aux = sqldf(query) #df
GSEIDs = df_aux$GSEID  # character

### Creamos la tabla GSE | GSM | Sample_Source
## Obtenemos GSMs y samples sources

for (i in 1:length(GSEIDs)) {
  print(i)
  GSE = GSEIDs[i]
  # extraemos la información disponible del estudio en GEOMetadb:
  sample_query <- paste("SELECT DISTINCT gse.gse, gsm.gsm, gsm.gpl,gsm.source_name_ch1 as sample_source, gsm.organism_ch1 as organism",
                        " FROM",
                        " gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm",
                        " JOIN gse ON gse_gsm.gse=gse.gse",
                        " WHERE",
                        " gse.gse ='",GSE,"'", sep ='')
  
  rs_sq <- dbGetQuery(con, sample_query)
  print(GSE)
  print(length(rs_sq$gsm))
  n = length(rs_sq$gsm)
  if (n == 0) {
    # Si no hay información de las muestras el estudio no está en la imagen de GEOmetadb
    # Estudio no detectado: Completamos manualmente
    gse = c(GSE)
    gsm = c("NA")
    gpl = c(GPLID)
    sample_source = c("NA")
    organism = c("NA")
    rs_sq = data.frame(gse,gsm,gpl,sample_source,organism)
  }
  if (i == 1) {
    gse_gsm_relational = rs_sq
    conteos = c(n)
}
  else{
    #Append rs_sq (salida de la query) to res_df (dataframe final), añadir debajo
    gse_gsm_relational = rbind(gse_gsm_relational, rs_sq)
    conteos = c(conteos, n)
  }
}
# Tiene que ser igual al nº de estudios:
unique(gse_gsm_relational$gse)

tabla_conteos_gse_gsm = data.frame(GSEIDs)
tabla_conteos_gse_gsm$sample_count = conteos

## guardamos la relación gse - nmuestras
write.csv2(tabla_conteos_gse_gsm,paste0('./',GPLID,'/',GPLID,'_gse_sampleCount.csv'), row.names = FALSE)

## guardamos la relación gse_gsm_relational
write.csv2(gse_gsm_relational,paste0('./',GPLID,'/',GPLID,'_gse_gsm_relational.csv'), row.names = FALSE)

#write.csv2(gse_gsm_relational,'gse_gsm_relational.csv', row.names = FALSE)
############### PARTE 2 ###########################
# PROCESAMOS gse_gsm_relational para identificar las muestras que se corresponden con tejido adiposo en cada estudio.
# Contamos numero de muestras totales
# Contamos numero de muestras de tejido adiposo - '%adipo%' xq puede ser adipose/adipocyte/...
# Si el numero de muestras no es el mismo

## Añadimos esta información a la tabla resumen de los adipose_GSEs
## GSE | GSM | SAMPLE_SOURCE | indicador_MTej
## SELECT GSE, GSM FROM table WHERE SAMPLE_SOURCE LIKE '%adipose%'

## necesito un indicador de si el estudio es monotejido o multitejido

## COUNT (GSM) FROM table WHERE gse = 'GSEXXXXX' ; numero total de muestras por estudio
## COUN (GSM) FROM table WHERE gse = 'GSEXXXXX' AND sample_source LIKE '%adipose%'; numero de muestras de tejido adiposo

### GSE | nsamples_total | nsamples_adipose -> si nsamples total - nsamples_adipose != 0 -> estudio multitejido
# Leemos el fichero gse_gsm_relational, con la relacion de estudios y muestras
GPLID = "GPL16570"
gse_gsm_relational <- read_delim(paste0(GPLID,"/",GPLID,"_gse_gsm_relational.csv"),
";", escape_double = FALSE, trim_ws = TRUE)

indicadorMT = c() # Vector para almacenar el indicador monotejido/multitejido
nSamples_tot = c()
nSamples_adi = c()
adiGSMs = c()
organism = "Mus musculus" # Homo sapiens

GSEIDs = unique(gse_gsm_relational$gse)
for (i in 1:length(GSEIDs)) {# GPL_Express/brain_GPLs
  print(i)
  GSE = GSEIDs[i]
  count_allGPLsamples = paste('SELECT COUNT(gsm) FROM gse_gsm_relational WHERE gse = "',GSE,'" AND',
                           ' gpl = "',GPLID,'" AND',
                           ' organism LIKE "%',organism,'%"', sep='')
  
  count_adiposeGPLsamples = paste('SELECT COUNT(gsm) FROM gse_gsm_relational WHERE gse = "',GSE,'" AND',
                               ' gpl = "',GPLID,'" AND',
                               ' sample_source LIKE "%adipo%" AND',
                               ' organism LIKE "%',organism,'%"', sep='')
  
  df_AllS = sqldf(count_allGPLsamples)
  df_AdiS = sqldf(count_adiposeGPLsamples)
  
  nTot = df_AllS$`COUNT(gsm)`
  nAdi = df_AdiS$`COUNT(gsm)`
  
  nSamples_tot = c(nSamples_tot, nTot)
  nSamples_adi = c(nSamples_adi, nAdi)
  
  if (df_AllS$`COUNT(gsm)` != df_AdiS$`COUNT(gsm)`) {
    ## Indicador de estudio multitejido:
    ind = 'MULTITEJIDO'
    indicadorMT = c(indicadorMT, ind)
    ## Queremos los identificadores de las muestras de cerebro
    get_adipose_samples = paste('SELECT gsm FROM gse_gsm_relational WHERE GSE = "',GSE,'" AND',
                                ' gpl = "',GPLID,'" AND',
                                ' sample_source LIKE "%adipo%" AND',
                                ' organism LIKE "%',organism,'%"', sep='')
    adisamples = sqldf(get_adipose_samples)
    
    # Extraemos los GSMs
    GSM_aux = adisamples$gsm
    # Convertimos la lista en una cadena:
    q = GSM_aux
    GSM_string = ''
    for (i in 1:length(q)) {
      gsm = q[i]
      if (i == 1) { # sabemos que GPL570 no es la primera
        gsm = paste("\'",gsm,"\'", sep ='')  # necesitamos añadirle las comillas
        GSM_string = paste(GSM_string, gsm, sep = '' )
      }
      else {
        gsm = paste("\'",gsm,"\'", sep ='')  # necesitamos añadirle las comillas
        GSM_string = paste(GSM_string, gsm, sep = ', ' )
      }
    }
    adiGSMs = c(adiGSMs, GSM_string)
  } else if(df_AllS$`COUNT(gsm)` == 0|| df_AdiS$`COUNT(gsm)` == 0){
    ind = "NO INFO"
    indicadorMT = c(indicadorMT, ind)
    adiGSMs = c(adiGSMs,'NULL' )  
    }else{
    ## Indicador de estudio monotejido
    ind = 'MONOTEJIDO'
    indicadorMT = c(indicadorMT, ind)
    ## No necesitamos los identificadores de las muestras
    adiGSMs = c(adiGSMs,'ALL' )
  }
}

summary = data.frame(GSEIDs,nSamples_tot,nSamples_adi,indicadorMT, adiGSMs)
write.csv(summary,file = paste0('./',GPLID,'/',GPLID,'_summary.csv'), row.names = FALSE)

#####################################################################################
#####################################################################################
## LAS MUESTRAS DE UN ESTUDIO NO TIENEN PORQUE SER LAS MISMAS PARA TODAS LAS
## PLATAFORMAS DEL ESTUDIO.
## Nº de muestras total del estudio
## Nº de muestras de la plataforma de interés en el estudio
## Nº de muestras de la plataforma de interés que sean de tej.adiposo en el estudio
#####################################################################################
#####################################################################################



