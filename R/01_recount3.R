# CÓDIGO PARA SELECCIONAR EL PROYECTO USANDO RECOUNT3
# Instalar paquete
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("recount3", "SummarizedExperiment", "GenomicRanges"))

# Cargar los paquetes
library(recount3)
library(SummarizedExperiment)
library(GenomicRanges)

# Obtener la lista de proyectos disponibles
human_projects <- available_projects()

# Ver los proyectos disponibles
dim(human_projects)

# Esto nos indica cuántos proyectos están disponibles (número de filas)
# y cuántas columnas de información se proporcionan para cada proyecto.

# Mostrar las primeras filas para inspeccionar su estructura y contenido
head(human_projects)

# Seleccionar un estudio de interés
human_projects[709, ]

# Filtrar el dataframe para seleccionar un proyecto específico basado en su ID y tipo
project_info <- subset(
  human_projects,
  project == "SRP075398" & project_type == "data_sources"
)

# Mostrar la información del proyecto seleccionado para confirmar que se ha
# filtrado correctamente
project_info

# Crear un objeto de tipo RangedSummarizedExperiment (RSE) con la información a nivel de genes
rse_gene_SRP075398 <- create_rse(project_info)

# Explorar el objeto RSE
rse_gene_SRP075398

## Información sobre el RSE creado
metadata(rse_gene_SRP075398)

## Número de genes y número de muestras
dim(rse_gene_SRP075398)

# Información sobre los genes
rowRanges(rse_gene_SRP075398)
