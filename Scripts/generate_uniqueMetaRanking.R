library(readr)
library(tibble)
distance_CV <- read_csv("Human/a/Proximidad/distance_CV.csv")
distance_IQRm <- read_csv("Human/a/Proximidad/distance_IQRmedian.csv")
distance_MADm <- read_csv("Human/a/Proximidad/distance_MADmedian.csv")
exprs = read.delim("Human/a/Proximidad/Unify/gtex_median_tpm.txt", header = TRUE)

# Ordenamos alfabeticamente:
distance_CV <- distance_CV[order(distance_CV$SYMBOL),]
distance_IQRm <- distance_IQRm[order(distance_IQRm$SYMBOL),]
distance_MADm <- distance_MADm[order(distance_MADm$SYMBOL),]

# Extraemos los genes:
genes_cv = distance_CV$SYMBOL
genes_iqr = distance_IQRm$SYMBOL
genes_mad = distance_MADm$SYMBOL

# Nos aseguramos que se corresponden:
table(genes_cv == genes_iqr)
table(genes_iqr == genes_mad)

#Extraemos los identificadores:
SYMBOL = genes_cv
GENENAME = distance_CV$GENENAME
stats = c('CV','IQRm','MADm')

# Extraemos las posiciones de estabilidad en hombres:  
rank_cv = distance_CV$men_Ranking
rank_iqr = distance_IQRm$men_Ranking
rank_mad = distance_MADm$men_Ranking

# Unificamos las posiciones de los estadisticos en un solo indicador promedio:
df_men = data.frame(SYMBOL,rank_cv, rank_iqr, rank_mad)  
df_men$rank_promedio = round(((df_men$rank_cv + df_men$rank_iqr + df_men$rank_mad)/3),2)
men_promedio = df_men$rank_promedio
#write.csv(df_men,paste0("./Unify/ranking_men.csv"), row.names = FALSE)


# Extraemos las posiciones de estabilidad en mujeres:
rank_cv = distance_CV$women_Ranking
rank_iqr = distance_IQRm$women_Ranking
rank_mad = distance_MADm$women_Ranking

# Unificamos las posiciones de los estadisticos en un solo indicador promedio:
df_women = data.frame(SYMBOL, rank_cv, rank_iqr, rank_mad)  
df_women$rank_promedio = round(((df_women$rank_cv + df_women$rank_iqr + df_women$rank_mad)/3),2)
women_promedio = df_women$rank_promedio
#write.csv(df_women,paste0("./Unify/ranking_women.csv"), row.names = FALSE)


# Creamos el df conjunto:
comparative = data.frame(SYMBOL, GENENAME ,men_promedio, women_promedio)
comparative$distancia = comparative$men_promedio - comparative$women_promedio

# Incorporamos la información de expresión con un join:
final_df = merge(x = comparative, y = exprs, by ="SYMBOL", all.x = TRUE)
length(unique(final_df$SYMBOL))
colnames(final_df)[2] <- "GENENAME"


# Reordenamos las columnas y eliminamos la información que no queremos:
final_df2 = final_df[,c(1,2,6,9,3,4,5,7)]


# Anotamos:
functional_annotation_CC <- read_csv("functional_annotation_CC.csv")
functional_annotation_MF <- read_csv("functional_annotation_MF.csv")
functional_annotation_BP <- read_csv("functional_annotation_BP.csv")

ONTOLOGIES = c("BP", "MF", "CC")

for (o in 1:length(ONTOLOGIES)){
  ONT = ONTOLOGIES[o]
  df_aux = final_df2
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

  write.csv(ndf,paste0("./Unify/comparative_annotated_",ONT,".csv"), row.names = FALSE)

  }



  
  