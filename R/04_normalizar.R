# CÓDIGO PARA NORMALIZAR LOS DATOS

# Crear un objeto DGEList para normalización
dge <- DGEList(
  counts = assay(rse_gene_SRP075398, "counts"),
  genes = rowData(rse_gene_SRP075398)
)

# Normalización TMM
dge <- calcNormFactors(dge)

dge

