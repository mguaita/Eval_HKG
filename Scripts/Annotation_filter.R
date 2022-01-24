library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(readr)
library(sqldf)

## Comprobamos la información de anotación:
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(GO.db)
help("SYMBOL")

### ANOTACION ###
# 1. Extraemos el listado de genes a anotar
RP_FEMALE_medianCV <- read_csv("CV/RP_FEMALE_medianCV.csv")
RP_MALE_medianCV <- read_csv("CV/RP_MALE_medianCV.csv")
RP_FEMALE_medianIQRmedian <- read_csv("IQR/RP_FEMALE_medianIQRmedian.csv")
RP_MALE_medianIQRmedian <- read_csv("IQR/RP_MALE_medianIQRmedian.csv")
RP_FEMALE_medianMADmedian <- read_csv("MAD/RP_FEMALE_medianMADmedian.csv")
RP_MALE_medianMADmedian <- read_csv("MAD/RP_MALE_medianMADmedian.csv")


#Cogemos un listado de referencia:
gene_list = as.vector(RP_FEMALE_medianCV$Gene) # 41975 genes

# Corregimos anotación RNA18S (estandarizamos la nom)
#gene_list = replace(gene_list, gene_list=="RNA18S", "RNA18SN5")

# 2. Convertimos: Gene Name es la descripción del gen
annotated <- select(org.Hs.eg.db, keys=gene_list, columns=c("GENENAME"), keytype="SYMBOL") #41979 !COMPROBAR RPS18#
table(is.na(annotated$GENENAME))
# FALSE  TRUE 
# 25214 16765 ; tenemos 16765 genes sin anotación
table(is.na(annotated$SYMBOL))

#annotated <- merge(annotated, gene_name, by.x = "SYMBOL", by.y = "Gene name", all.x = TRUE)

## Tenemos duplicados:
annotated[duplicated(annotated$SYMBOL),] 
#SYMBOL                                  GENENAME
#17801  MEMO1               mediator of cell motility 1
#31846    TEC transient erythroblastopenia of childhood
#32173   MMD2              Miyoshi muscular dystrophy 2
#39758    HBD             hypophosphatemic bone disease

## Los limpiamos manualmente, son solo 4 casos:
annotated[which(annotated$SYMBOL == "MEMO1"),] # Nos quedamos con la 2da acepción 17801, eliminamos fila 17800
annotated[which(annotated$SYMBOL == "TEC"),] # eliminamos 31846
annotated[which(annotated$SYMBOL == "MMD2"),] # eliminamos 32173
annotated[which(annotated$SYMBOL == "HBD"),] # eliminamos 39758

## Eliminamos las filas que no nos interesan:
annotated = annotated[-c(17800,31846,32173,39758),] # 41975


### Anotamos con el fichero de información de Biomart (creo que es más completa):
# Gene description, gene name, y terminos GO (ID, term name y term definition)


### FILTRO DE GENES ###
# Excluimos los genes que contengan en su descripción
clean_data <- sqldf("Select * FROM annotated
where GENENAME NOT Like '%pseudogene%' 
AND GENENAME NOT Like '%microRNA%' 
AND GENENAME NOT Like '%small nuclear%' 
AND GENENAME NOT Like '%small nucleolar%' 
AND GENENAME NOT Like '%non-protein coding%'
AND GENENAME NOT Like '%long non-coding%'
AND GENENAME NOT Like '%uncharacterized LOC%'")

table(is.na(clean_data)) 
# Listado de genes limpio:
annotated_genes = as.vector(clean_data$SYMBOL) # 18973

### RESULTADO: 
## 41975 genes iniciales
## 25214 genes anotados con AnnotationDbi
## 18973 genes superan el filtro


### ANOTACION FUNCIONAL ### 
#a = "RNA18S5"
# 1 Extraemos la información de los genes limpios - 18973
## Tenemos 2 dfs, uno con los genes limpios y los identificadores ID GO asociados, y
## otro con la información de los terminos GO (definicion y nombre)
GO_df <- select(org.Hs.eg.db, keys=annotated_genes, columns=c("GO"), keytype="SYMBOL") #316940
genes_with_GO_associated = unique(GO_df$SYMBOL) # tenemos información de todos los genes 18973

table(is.na(GO_df$GO)) # Tenemos nulos
GO_df = GO_df[!is.na(GO_df$GO),]

GO_IDs <- GO_df$GO
GO_info <- select(GO.db, keys=GO_IDs, columns=c("DEFINITION", "TERM", "ONTOLOGY"), keytype="GOID") #316939
table(is.na(GO_info$DEFINITION))
GO_info = GO_info[!is.na(GO_info$TERM),]


# Procesamos cada ontologia por separado:
BP_df = GO_df[which(GO_df$ONTOLOGY == "BP"),] #148500
MF_df = GO_df[which(GO_df$ONTOLOGY == "MF"),] #75832
CC_df = GO_df[which(GO_df$ONTOLOGY == "CC"),] #91295

BP_info = GO_info[which(GO_info$ONTOLOGY == "BP"),] #148500
MF_info = GO_info[which(GO_info$ONTOLOGY == "MF"),] #75832
CC_info = GO_info[which(GO_info$ONTOLOGY == "CC"),] #91295


#setwd("/clinicfs/userhomes/mguaita/Human/adipose/Automated/RESULTS/Proximidad")
## Anotamos los GO IDs (GO_df) de los genes con su correspondiente información (GO_info):
func_annot_BP <- sqldf("Select DISTINCT SYMBOL, GO, BP_info.TERM, BP_info.DEFINITION, BP_info.ONTOLOGY FROM BP_df
                     JOIN BP_info ON BP_df.GO == BP_info.GOID")
write.csv(func_annot_BP,paste0("functional_annotation_BP.csv"), row.names = FALSE)

func_annot_MF <- sqldf("Select DISTINCT SYMBOL, GO, MF_info.DEFINITION, MF_info.TERM, MF_info.ONTOLOGY FROM MF_df
                     JOIN MF_info ON MF_df.GO == MF_info.GOID")
write.csv(func_annot_MF,paste0("functional_annotation_MF.csv"), row.names = FALSE)

func_annot_CC <- sqldf("Select DISTINCT SYMBOL, GO, CC_info.DEFINITION, CC_info.TERM, CC_info.ONTOLOGY FROM CC_df
                     JOIN CC_info ON CC_df.GO == CC_info.GOID")
write.csv(func_annot_CC,paste0("functional_annotation_CC.csv"), row.names = FALSE)


## Ya tenemos los genes anotados:
write.csv(func_annot_BP,paste0("functional_annotation_BP.csv"), row.names = FALSE)
write.csv(func_annot_MF,paste0("functional_annotation_MF.csv"), row.names = FALSE)
write.csv(func_annot_CC,paste0("functional_annotation_CC.csv"), row.names = FALSE)

### CALCULO INDICADOR DE PROXIMIDAD ###
## Estandarizamos la nomenclatura:
# RP_FEMALE_medianCV$Gene[RP_FEMALE_medianCV$Gene == "RNA18S"] <- "RNA18SN5"
# RP_MALE_medianCV$Gene[RP_MALE_medianCV$Gene == "RNA18S"] <- "RNA18SN5"
# 
# RP_FEMALE_medianIQRmedian$Gene[RP_FEMALE_medianIQRmedian$Gene == "RNA18S"] <- "RNA18SN5"
# RP_MALE_medianIQRmedian$Gene[RP_MALE_medianIQRmedian$Gene == "RNA18S"] <- "RNA18SN5"
# 
# RP_FEMALE_medianMADmedian$Gene[RP_FEMALE_medianMADmedian$Gene == "RNA18S"] <- "RNA18SN5"
# RP_MALE_medianMADmedian$Gene[RP_MALE_medianMADmedian$Gene == "RNA18S"] <- "RNA18SN5"



# 1- Extraemos los genes anotados y sus valores de RP (H/M)
# 2- Ordenamos estos valores y generamos un nuevo metaRanking (H/M)
# 3- Calculamos el indicador de proximidad en base a la posición que ocupan en el ranking (H - M)

## 1 Extraemos genes y valores RP del listado de genes anotados limpios (H/M) ##
stats = c("CV", "IQRmedian", "MADmedian")
for (i in 1: length(stats)){
  stat = stats[i]
  
  H_genes <- sqldf(paste0("Select SYMBOL, GENENAME, male.RP FROM clean_data
  JOIN RP_MALE_median",stat," AS male ON clean_data.SYMBOL = male.Gene"))
  names(H_genes)[3] <- "men_RP"
  
  M_genes <- sqldf(paste0("Select SYMBOL, GENENAME, female.RP FROM clean_data
  JOIN RP_FEMALE_median",stat," AS female ON clean_data.SYMBOL = female.Gene"))
  names(M_genes)[3] <- "women_RP"
  
  ## 2 Ordenamos de menor a mayor segun valor RP (H/M) ##
  H_genes = H_genes[order(H_genes$men_RP),]
  H_genes$men_Ranking <- seq.int(nrow(H_genes))
  
  M_genes = M_genes[order(M_genes$women_RP),]
  M_genes$women_Ranking <- seq.int(nrow(M_genes))
  
  
  ## 3 Unimos los rankings de hombres y mujeres y calculamos el estadistico de distancia
  comparative  <- sqldf("Select male.SYMBOL, male.GENENAME, male.men_Ranking, female.women_Ranking FROM H_genes as male
  JOIN M_genes AS female ON male.SYMBOL = female.SYMBOL")
  
  comparative$distance <- comparative$men_Ranking - comparative$women_Ranking
  comparative_2 = comparative[order(comparative$distance),]
  write.csv(comparative_2,paste0("distance_",stat,".csv"), row.names = FALSE)
}
## Interpretacion:
## diseño: Hombres - Mujeres
## Valores negativos hombre muy estable - mujeres variable
## Valores positivos mujeres muy estable - hombres variable
## Valores extremos indican niveles de estabilidad diferencialess
## Valores próximos a 0 indican niveles de estabilidad similares



# genes con mayor distancia en hombres (primeros 200)
h_CV = distance_CV[1:200,]
write.csv(h_CV,paste0("h_CV.csv"), row.names = FALSE)
# genes con mayor distancia en mujeres (ultimos 200)
m_CV = distance_CV[18773:18972,]
write.csv(m_CV,paste0("m_CV.csv"), row.names = FALSE)
# genes con estabilidad relativa similar (+-100 del 0)
n_CV = distance_CV[9338:9538,]
write.csv(n_CV,paste0("n_CV.csv"), row.names = FALSE)



# genes con mayor distancia en hombres (primeros 200)
h_IQR = distance_IQRmedian[1:200,] 
write.csv(h_IQR,paste0("h_IQR.csv"), row.names = FALSE)
# genes con mayor distancia en mujeres (ultimos 200)
m_IQR = distance_IQRmedian[18773:18972,]
write.csv(m_IQR,paste0("m_IQR.csv"), row.names = FALSE)
# genes con estabilidad relativa similar (+-100 del 0)
n_IQR = distance_IQRmedian[9185:9385,]
write.csv(n_IQR,paste0("m_IQR.csv"), row.names = FALSE)


# genes con mayor distancia en hombres (primeros 200)
h_MAD = distance_MADmedian[1:200,] 
write.csv(h_MAD,paste0("h_MAD.csv"), row.names = FALSE)
# genes con mayor distancia en mujeres (ultimos 200)
m_MAD = distance_MADmedian[18773:18972,]
write.csv(m_MAD,paste0("m_MAD.csv"), row.names = FALSE)
# genes con estabilidad relativa similar (+-100 del 0)
n_MAD = distance_MADmedian[9391:9591,]
write.csv(n_MAD,paste0("m_MAD.csv"), row.names = FALSE)



### ANOTAR RESULTADOS FINALES ###
functional_annotation_CC <- read_csv("functional_annotation_CC.csv")
functional_annotation_MF <- read_csv("functional_annotation_MF.csv")
functional_annotation_BP <- read_csv("functional_annotation_BP.csv")

distance_CV <- read_csv("distance_CV.csv")
distance_MADmedian <- read_csv("distance_MADmedian.csv")
distance_IQRmedian <- read_csv("distance_IQRmedian.csv")

stats = c("CV", "IQRmedian", "MADmedian")
ONTOLOGIES = c("BP", "MF", "CC")
for ( s in 1:length(stats)){
  stat = stats[s]
  for (o in 1:length(ONTOLOGIES)){
    ONT = ONTOLOGIES[o]
    df_aux = get(paste0("distance_",stat))
    gene_list = as.vector(df_aux$SYMBOL)
    fc = get(paste0("functional_annotation_",ONT))
    all_GO_terms = c()
    all_GO_IDs = c()
    for (i in 1:length(gene_list)){
      gene = gene_list[i] # character
      df_aux2 = fc[which(fc$SYMBOL == gene),]  # df
      GO_terms_v = df_aux2$TERM # vector
      GO_terms_s = paste(GO_terms_v, collapse = ", ")
      all_GO_terms = c(all_GO_terms,GO_terms_s)
      GO_IDs_v = df_aux2$GO # vector
      GO_IDs_s = paste(GO_IDs_v, collapse = ", ")
      all_GO_IDs = c(all_GO_IDs, GO_IDs_s)
      }
    
    ndf = df_aux
    ndf$GO_terms = all_GO_terms
    ndf$GO_IDs = all_GO_IDs
    if (stat == "IQRmedian"){
      stat = "IQRm"
    }
    if (stat == "MADmedian"){
      stat = "MADm"
    }
    write.csv(ndf,paste0("./fully_annotated/distance_",stat,"_annot",ONT,".csv"), row.names = FALSE)
    if (stat == "IQRm"){
      stat = "IQRmedian"
    }
    if (stat == "MADm"){
      stat = "MADmedian"
    }
}
}


































folder = getwd()    # path to folder that holds multiple .csv files
file_list <- list.files(path=paste0(folder,"/"), pattern="*.csv") # create list of all .csv files in folder

# read in each .csv file in file_list and create a data frame with the same name as the .csv file
for (i in 1:length(file_list)){
  assign(file_list[i], 
         read.csv(paste(folder,"/", file_list[i], sep=''))
  )}

CV_BP = distance_CV_annotBP.csv[9138:9738,]
write.csv(CV_BP,paste0("./similar_values/comunes_CV_BP.csv"), row.names = FALSE)

CV_CC = distance_CV_annotCC.csv[9138:9738,]
write.csv(CV_CC,paste0("./similar_values/comunes_CV_CC.csv"), row.names = FALSE)

CV_MF = distance_CV_annotMF.csv[9138:9738,]
write.csv(CV_MF,paste0("./similar_values/comunes_CV_MF.csv"), row.names = FALSE)


IQR_BP = distance_IQRm_annotBP.csv[8985:9585,]
write.csv(IQR_BP,paste0("./similar_values/comunes_IQR_BP.csv"), row.names = FALSE)

IQR_CC =distance_IQRm_annotCC.csv[8985:9585,]
write.csv(IQR_CC,paste0("./similar_values/comunes_IQR_CC.csv"), row.names = FALSE)

IQR_MF = distance_IQRm_annotMF.csv[8985:9585,]
write.csv(IQR_MF,paste0("./similar_values/comunes_IQR_MF.csv"), row.names = FALSE)

MAD_BP = distance_MADm_annotBP.csv[9191:9791,]
write.csv(MAD_BP,paste0("./similar_values/comunes_MAD_BP.csv"), row.names = FALSE)

MAD_CC = distance_MADm_annotCC.csv[9191:9791,]
write.csv(MAD_CC,paste0("./similar_values/comunes_MAD_CC.csv"), row.names = FALSE)

MAD_MF = distance_MADm_annotMF.csv[9191:9791,]
write.csv(MAD_MF,paste0("./similar_values/comunes_MAD_MF.csv"), row.names = FALSE)










BiocManager::install("ReportingTools")
library(ReportingTools)

# Nombre del output:
edgeR_output = "edgeR.DE"

# Generación del informe
htmlEDGERreport = HTMLReport(shortName = edgeR_output, title = "edgeR_DE",
                             reportDirectory = "./reports")

# Acabamos todo
publish(def_report_edgeR,htmlEDGERreport)
finish(htmlEDGERreport)









#### COMMON MOST STABLE GENES ###

distance_CV <- read_csv("distance_CV.csv")
distance_IQRmedian <- read_csv("distance_IQRmedian.csv")
distance_MADmedian <- read_csv("distance_MADmedian.csv")

# Prueba con CV: (Repetimos codigo con los otros 2 stats)
# 1.- Ordenamos los genes por su estabilidad en hombres - extraemos los 900 mas estables (segun el orden)
datos <- distance_CV[with(distance_CV, order(distance_CV$men_Ranking)), ]

#  Extraemos los genes más estables de cada sexo:
men_NH = as.vector(datos$SYMBOL) # 41975 genes
men_NH = men_NH[1:900]


datos <- distance_CV[with(distance_CV, order(distance_CV$women_Ranking)), ]

women_NH = as.vector(datos$SYMBOL) # 41975 genes
women_NH = women_NH[1:900]

# 2.- Buscamos los comunes:
length(intersect(men_NH, women_NH))
common_genes = intersect(men_NH, women_NH) # 203 genes comunes

# 3.- Extraemos sus posiciones del DF de las distancias:
common_genes_string = ''
for (i in 1:length(common_genes)) {
  g = common_genes[i]
  if (i == 1) { # sabemos que GPL570 no es la primera
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = '' )
  }
  else {
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = ', ' )
  }
}


common_df = sqldf(paste0("select SYMBOL, GENENAME, men_Ranking, women_Ranking, distance ",
                         "FROM distance_CV ",
                         "WHERE SYMBOL IN (",common_genes_string,")"))

write.csv(common_df,paste0("./Common/common_CV.csv"), row.names = FALSE)



#### FUNCION DE ANOTACION ####


### ANOTAR RESULTADOS FINALES ###
functional_annotation_CC <- read_csv("functional_annotation_CC.csv")
functional_annotation_MF <- read_csv("functional_annotation_MF.csv")
functional_annotation_BP <- read_csv("functional_annotation_BP.csv")

common_CV <- read_csv("./Common/common_CV.csv") # 203 genes
common_MADmedian <- read_csv("./Common/common_MADmedian.csv") # 224 genes
common_IQRmedian <- read_csv("./Common/common_IQRmedian.csv") # 172 genes

stats = c("CV", "IQRmedian", "MADmedian")
ONTOLOGIES = c("BP", "MF", "CC")
for ( s in 1:length(stats)){
  stat = stats[s]
  for (o in 1:length(ONTOLOGIES)){
    ONT = ONTOLOGIES[o]
    df_aux = get(paste0("common_",stat))
    gene_list = as.vector(df_aux$SYMBOL)
    fc = get(paste0("functional_annotation_",ONT))
    all_GO_terms = c()
    all_GO_IDs = c()
    for (i in 1:length(gene_list)){
      gene = gene_list[i] # character
      df_aux2 = fc[which(fc$SYMBOL == gene),]  # df
      GO_terms_v = df_aux2$TERM # vector
      GO_terms_s = paste(GO_terms_v, collapse = ", ")
      all_GO_terms = c(all_GO_terms,GO_terms_s)
      GO_IDs_v = df_aux2$GO # vector
      GO_IDs_s = paste(GO_IDs_v, collapse = ", ")
      all_GO_IDs = c(all_GO_IDs, GO_IDs_s)
    }
    
    ndf = df_aux
    ndf$GO_terms = all_GO_terms
    ndf$GO_IDs = all_GO_IDs
    if (stat == "IQRmedian"){
      stat = "IQRm"
    }
    if (stat == "MADmedian"){
      stat = "MADm"
    }
    write.csv(ndf,paste0("./Common/common",stat,"_annot",ONT,".csv"), row.names = FALSE)
    if (stat == "IQRm"){
      stat = "IQRmedian"
    }
    if (stat == "MADm"){
      stat = "MADmedian"
    }
  }
}







###### CODIGO REPETIDO ####
# 1.- Ordenamos los genes por su estabilidad en hombres - extraemos los 900 mas estables (segun el orden)
datos <- distance_IQRmedian[with(distance_IQRmedian, order(distance_IQRmedian$men_Ranking)), ]

#  Extraemos los genes más estables de cada sexo:
men_NH = as.vector(datos$SYMBOL) # 41975 genes
men_NH = men_NH[1:900]


datos <- distance_IQRmedian[with(distance_IQRmedian, order(distance_IQRmedian$women_Ranking)), ]

women_NH = as.vector(datos$SYMBOL) # 41975 genes
women_NH = women_NH[1:900]

# 2.- Buscamos los comunes:
length(intersect(men_NH, women_NH))
common_genes = intersect(men_NH, women_NH) # 203 genes comunes

# 3.- Extraemos sus posiciones del DF de las distancias:
common_genes_string = ''
for (i in 1:length(common_genes)) {
  g = common_genes[i]
  if (i == 1) { # sabemos que GPL570 no es la primera
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = '' )
  }
  else {
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = ', ' )
  }
}


common_df = sqldf(paste0("select SYMBOL, GENENAME, men_Ranking, women_Ranking, distance ",
                         "FROM distance_IQRmedian ",
                         "WHERE SYMBOL IN (",common_genes_string,")"))

write.csv(common_df,paste0("./Common/common_IQRmedian.csv"), row.names = FALSE)




########
# 1.- Ordenamos los genes por su estabilidad en hombres - extraemos los 900 mas estables (segun el orden)
datos <- distance_MADmedian[with(distance_MADmedian, order(distance_MADmedian$men_Ranking)), ]

#  Extraemos los genes más estables de cada sexo:
men_NH = as.vector(datos$SYMBOL) # 41975 genes
men_NH = men_NH[1:900]


datos <- distance_MADmedian[with(distance_MADmedian, order(distance_MADmedian$women_Ranking)), ]

women_NH = as.vector(datos$SYMBOL) # 41975 genes
women_NH = women_NH[1:900]

# 2.- Buscamos los comunes:
length(intersect(men_NH, women_NH))
common_genes = intersect(men_NH, women_NH) # 172 genes comunes

# 3.- Extraemos sus posiciones del DF de las distancias:
common_genes_string = ''
for (i in 1:length(common_genes)) {
  g = common_genes[i]
  if (i == 1) { # sabemos que GPL570 no es la primera
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = '' )
  }
  else {
    g = paste("\'",g,"\'", sep ='')  # necesitamos añadirle las comillas
    common_genes_string = paste(common_genes_string, g, sep = ', ' )
  }
}


common_df = sqldf(paste0("select SYMBOL, GENENAME, men_Ranking, women_Ranking, distance ",
                         "FROM distance_MADmedian ",
                         "WHERE SYMBOL IN (",common_genes_string,")"))

write.csv(common_df,paste0("./Common/common_MADmedian.csv"), row.names = FALSE)




