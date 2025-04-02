# Script: analisis_st003741.R
# Descripción: funciones utilizadas para el análisis de datos metabolómicos del estudio ST003741

# -----------------------
# Cargar librerías necesarias
# -----------------------
library(readr)
library(dplyr)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(stringr)
library(ggplot2)
library(reshape2)
library(SummarizedExperiment)

# -----------------------
# Función para leer los datos .txt desde Metabolomics Workbench
# -----------------------
leer_datos <- function(nombre_archivo, ruta_base) {
  archivo <- read_lines(paste0(ruta_base, nombre_archivo))
  tabla_inicio <- which(stringr::str_count(archivo, "	") > 5)[1]
  datos_brutos <- archivo[tabla_inicio:length(archivo)]
  tabla <- read_tsv(paste(datos_brutos, collapse = "
"), show_col_types = FALSE)
  if ("Factors" %in% colnames(tabla)) {
    tabla <- tabla %>% select(-Factors)
  }
  tabla_rownames <- make.unique(tabla[[1]])
  rownames(tabla) <- tabla_rownames
  tabla <- tabla[,-1]
  matriz <- t(apply(tabla, 2, as.numeric))
  colnames(matriz) <- rownames(tabla)
  rownames(matriz) <- colnames(tabla)
  return(matriz)
}

# -----------------------
# Función para extraer metadatos desde los nombres de muestra
# -----------------------
extraer_metadata <- function(matriz) {
  muestras <- rownames(matriz)
  condiciones <- str_match(muestras, "(iRBC|RBC).*?(BR|DMSO)")
  df <- data.frame(Sample = muestras,
                   Infection = ifelse(condiciones[,2] == "iRBC", "Infected", "Non-infected"),
                   Treatment = ifelse(condiciones[,3] == "BR", "Bilirubin", "DMSO"),
                   stringsAsFactors = TRUE)
  rownames(df) <- muestras
  df <- df[complete.cases(df), ]
  return(df)
}

# -----------------------
# Función para graficar un PCA usando SummarizedExperiment
# -----------------------
graficar_pca_se <- function(se, titulo) {
  matriz <- t(assay(se))
  metadata <- as.data.frame(colData(se))
  matriz <- matriz[, colSums(is.na(matriz)) == 0]
  res_pca <- PCA(matriz, graph = FALSE)
  fviz_pca_ind(res_pca,
               habillage = metadata$Treatment,
               addEllipses = TRUE,
               palette = "jco",
               title = titulo)
}

# -----------------------
# Función para graficar un heatmap con los 5 metabolitos más variables
# -----------------------
graficar_heatmap_se <- function(se, titulo) {
  matriz <- t(assay(se))
  metadata <- as.data.frame(colData(se))
  varianzas <- apply(matriz, 2, var, na.rm = TRUE)
  top_vars <- names(sort(varianzas, decreasing = TRUE))[1:5]
  matriz_top <- matriz[, top_vars]
  matriz_top <- na.omit(matriz_top)
  matriz_top_scaled <- t(scale(t(matriz_top)))
  anotacion <- data.frame(Treatment = metadata$Treatment,
                          Infection = metadata$Infection)
  rownames(anotacion) <- rownames(matriz_top)
  pheatmap(matriz_top_scaled,
           annotation_row = anotacion,
           cluster_cols = TRUE,
           main = titulo)
}

# -----------------------
# Función para análisis exploratorio por grupos
# -----------------------
analisis_exploratorio_se <- function(se, titulo) {
  datos <- t(assay(se))
  metadata <- as.data.frame(colData(se))
  df <- as.data.frame(datos)
  df$Group <- interaction(metadata$Infection, metadata$Treatment)
  df_melt <- reshape2::melt(df, id.vars = "Group")
  ggplot(df_melt, aes(x = Group, y = value, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~variable, scales = "free_y") +
    theme_minimal() +
    labs(title = paste("Análisis exploratorio de metabolitos -", titulo), x = "Grupo", y = "Intensidad relativa")
}
