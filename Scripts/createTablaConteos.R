## CONEXION GEO
#Nos descargamos el esquema de la base de datos de GEO y establecemos la conexión:
library (GEOmetadb); packageDescription ("GEOmetadb", fields = "Version") #"1.44.0"
help (package = GEOmetadb)
# Cargamos el resto de paquetes:
library(sqldf)

###Get the database file: "GEOmetadb.sqlite"
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
getSQLiteFile()
file.info('GEOmetadb.sqlite')
dir ()


###Connect to the database - se crea la conexión (objeto conexión)
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
```


# BUCLE PARA GENERAR LA TABLA DE CONTEOS:
# Entrada: listado de plataformas
# Proceso: Va a la tabla, selecciona los estudios y los cuenta
# Salida: Para cada plataforma. número de estyuidos y sus identificadores

# 1) Cargamos el fichero gse_gpl_relational
library(readr)
gse_gpl_relational <- read_delim("gse_gpl_relational.csv", 
                                   +     "\t", escape_double = FALSE, trim_ws = TRUE)

# 2) Extraemos las plataformas unicas:

for (i in 1:length(gpls)) {# GPL_Express/brain_GPLs
  print(i)
  platform = gpls[i]
  
  # Construimos la query:
  query = paste('SELECT GSEID FROM gse_gpl_relational WHERE GPL = "',platform,'"', sep='')
  
  # Search:
  df_aux = sqldf(query) #df
  gseids = df_aux$GSEID  # character
  n = length(gseids)
  if (i ==1) {
    gse = list(gseids)
    conteos = c(n)}
  else{
    gse = append(gse, list(gseids))
    conteos = c(conteos, n)}
}

tabla_conteos = data.frame(gpls)
names(tabla_conteos)[names(tabla_conteos) == "gpls"] <- "gpl"
tabla_conteos$Nstudios = conteos
tabla_conteos$GSEIDs<- sapply(gse, paste0, collapse=",") 

write.csv2(tabla_conteos,'./tc_all_adipose.csv', row.names = FALSE)


write.csv2(full_brain,'./TABLACONTEOS_FULLBRAIN.csv', row.names = FALSE)


join <- merge(tabla_conteos, rs2, by.x=c('gpls'), by.y=c('gpl'), all= TRUE) 

### OBTENEMOS LA DESCRIPCIÓN DE LAS PLATAFORMAS
q = gpls
GPL_string = ''
for (i in 1:length(q)) {
  gpl = q[i]
  if (i == 1) { # sabemos que GPL570 no es la primera
    gpl = paste("\'",gpl,"\'", sep ='')  # necesitamos añadirle las comillas
    GPL_string = paste(GPL_string, gpl, sep = '' )
  }
  else {
    gpl = paste("\'",gpl,"\'", sep ='')  # necesitamos añadirle las comillas
    GPL_string = paste(GPL_string, gpl, sep = ', ' )
  }
}

mq2 <- paste("SELECT gpl.gpl, gpl.title",
  " FROM",
  " gpl",
  " WHERE",
  " gpl.gpl IN (",GPL_string,")", sep ="")
rs2 <- dbGetQuery(con, mq2)
View(rs2)

### unimos los df; ordenamos las columnas, escribimos el fichero
join <- merge(tabla_conteos, rs2, by.x=c('gpl'), by.y=c('gpl'), all= TRUE) 
tc2 = join[,c(1,4,2,3)]
write.csv2(tc2,'./adiposeGPLAnnotation.csv', row.names = FALSE)

