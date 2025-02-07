# CÓDIGO PARA PREPARAR LOS DATOS

# Convertir las cuentas por nucleotido a cuentas por lectura usando compute_read_counts()
assay(rse_gene_SRP075398, "counts") <- compute_read_counts(rse_gene_SRP075398)

# Inspeccionar la información experimental de cada muestra
rse_gene_SRP075398$sra.sample_attributes[]

# Expandir los atributos en columnas separadas para facilitar su uso
rse_gene_SRP075398 <- expand_sra_attributes(rse_gene_SRP075398)

# Extraer y mostrar las columnas que contienen atributos
colData(rse_gene_SRP075398)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP075398)))
]

# Ajustar el tipo de dato de las variables categóricas

rse_gene_SRP075398$sra_attribute.cell_line <- factor(rse_gene_SRP075398$sra_attribute.cell_line)

rse_gene_SRP075398$sra_attribute.source_name <- factor(tolower(rse_gene_SRP075398$sra_attribute.source_name))

rse_gene_SRP075398$sra_attribute.transfection <- factor(rse_gene_SRP075398$sra_attribute.transfection)

# Resumen estadístico de las variables seleccionadas
summary(as.data.frame(colData(rse_gene_SRP075398)[
  ,
  grepl("^sra_attribute.[cell_line|source_name|transfection]", colnames(colData(rse_gene_SRP075398)))
]))

# Calcular la proporción de lecturas asignadas a genes para evaluar la calidad de las muestras
rse_gene_SRP075398$assigned_gene_prop <-
  rse_gene_SRP075398$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP075398$recount_qc.gene_fc_count_all.total

# Resumen de la nueva variable para identificar si las muestras tienen una asignación adecuada de lecturas (valores cercanos a 1 son indicadores de buena calidad)
summary(rse_gene_SRP075398$assigned_gene_prop)
