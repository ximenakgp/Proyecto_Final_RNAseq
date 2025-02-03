# Cargar el paquete de R que incluye a SummarizedExperiment y todas las demás dependencias
library("recount3")

# Identificar el proyecto con el que le interese trabajar
human_projects <- available_projects()
# 2025-02-02 20:29:40.138416 caching file sra.recount_project.MD.gz.
# 2025-02-02 20:29:40.915193 caching file gtex.recount_project.MD.gz.
# 2025-02-02 20:29:41.625141 caching file tcga.recount_project.MD.gz.
dim(human_projects)
# [1] 8742    6
head(human_projects)
#   project organism file_source     project_home project_type n_samples
# 1 SRP107565    human         sra data_sources/sra data_sources       216
# 2 SRP149665    human         sra data_sources/sra data_sources         4
# 3 SRP017465    human         sra data_sources/sra data_sources        23
# 4 SRP119165    human         sra data_sources/sra data_sources         6
# 5 SRP133965    human         sra data_sources/sra data_sources        12
# 6 SRP096765    human         sra data_sources/sra data_sources         7

# Seleccionar un estudio de interés
human_projects[57, ]
#    project organism file_source     project_home project_type n_samples
# 57 SRP068565    human         sra data_sources/sra data_sources        20

## Colocar el ID del proyecto
project_info <- subset(
  human_projects,
  project == "SRP068565" & project_type == "data_sources"
)
project_info

#    project organism file_source     project_home project_type n_samples
# 57 SRP068565    human         sra data_sources/sra data_sources        20

# Crear un objeto de tipo RangedSummarizedExperiment (RSE) con la información a nivel de genes
rse_gene_SRP068565 <- create_rse(project_info)

## create_rse()es una función para GENCODE v26 (la anotación predeterminada para archivos humanos)

# Explorar el objeto RSE
rse_gene_SRP068565

# class: RangedSummarizedExperiment
# dim: 63856 20
# metadata(8): time_created recount3_version ... annotation recount3_url
# assays(1): raw_counts
# rownames(63856): ENSG00000278704.1 ENSG00000277400.1 ...
# ENSG00000182484.15_PAR_Y ENSG00000227159.8_PAR_Y
# rowData names(10): source type ... havana_gene tag
# colnames(20): SRR3105695 SRR3105677 ... SRR3105694 SRR3105696
# colData names(175): rail_id external_id ...
# recount_pred.curated.cell_line BigWigURL

## Información sobre el RSE creado
metadata(rse_gene_SRP068565)

# El estudio SRP068565 se compuso de 20 muestras, para las cuales tenemos 63,856 genes en GENCODE v26.
# La información específica de la anotación está disponible rowRanges()como se muestra a continuación con
# la columna gene_id utilizada para identificar genes en cada una de las anotaciones.

## Número de genes y número de muestras
dim(rse_gene_SRP068565)
#> [1] 63856    20

## Información sobre los genes
rowRanges(rse_gene_SRP068565)

# Preparar los datos para otras herramientas de análisis

# Convertir las cuentas por nucleotido a cuentas por lectura usando compute_read_counts().
assay(rse_gene_SRP068565, "counts") <- compute_read_counts(rse_gene_SRP068565)

# Hacer más fácil de usar la información del experimento
rse_gene_SRP068565 <- expand_sra_attributes(rse_gene_SRP068565)

colData(rse_gene_SRP068565)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP068565)))
]

# DataFrame with 20 rows and 2 columns
# sra_attribute.drug_treatment sra_attribute.source_name
# <character>               <character>
#   SRR3105695                         none                HeLa cells
# SRR3105677                    10 uM 1C8                HeLa cells
# SRR3105678                    10 uM 1C8                HeLa cells
# SRR3105679                    10 uM 1C8                HeLa cells
# SRR3105680                    10 uM 1C8                HeLa cells
# ...                                 ...                       ...
# SRR3105691                     5 uM 1C8                HeLa cells
# SRR3105692                     5 uM 1C8                HeLa cells
# SRR3105693                         none                HeLa cells
# SRR3105694                         none                HeLa cells
# SRR3105696                         none                HeLa cells


