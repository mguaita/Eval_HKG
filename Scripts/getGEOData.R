## FUNCION get_GEOData 
## Función de descarga de los datos de todos los estudios recogidos
## en una lista de identificadores GSE

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

