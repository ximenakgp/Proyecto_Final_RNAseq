# CÓDIGO PARA FILTRAR GENES DE BAJA EXPRESIÓN

library(edgeR)

# Guardar el objeto original
rse_gene_SRP075398_unfiltered <- rse_gene_SRP075398

# Visualizar la distribución de la proporción de lecturas asignadas a genes en cada muestra

hist(rse_gene_SRP075398$assigned_gene_prop,
     main = "Proporción de lecturas asignadas a genes",
     xlab = "Proporción asignada", col = "lightblue")

# Verificar si existen muestras de baja calidad antes del filtrado
table(rse_gene_SRP075398$assigned_gene_prop < 0.3)

# Filtrar las muestras con proporción de lecturas asignadas superior a 0.3
rse_gene_SRP075398 <- rse_gene_SRP075398[, rse_gene_SRP075398$assigned_gene_prop > 0.3]

# Crear un objeto DGEList, para el análisis diferencial usando edgeR
dge <- DGEList(counts = assay(rse_gene_SRP075398, "counts"))

# Filtrar genes de baja expresión considerando combinaciones de transfección y línea celular
keep <- filterByExpr(dge, group = interaction(
  rse_gene_SRP075398$sra_attribute.transfection,
  rse_gene_SRP075398$sra_attribute.cell_line
))
rse_gene_SRP075398 <- rse_gene_SRP075398[keep, ]

# Dimensiones finales
dim(rse_gene_SRP075398)

# Porcentaje de genes retenidos
round(nrow(rse_gene_SRP075398) / nrow(rse_gene_SRP075398_unfiltered) * 100, 2)
