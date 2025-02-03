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

