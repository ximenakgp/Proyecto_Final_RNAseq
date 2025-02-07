# CÓDIGO PARA REALIZAR EL ANÁLISIS DE EXPRESIÓN DIFERENCIAL

library(limma)
library(pheatmap)
library(RColorBrewer)

# Convertir los datos de conteo a valores log2 y ajusta las varianzas para hacerlos aptos para un análisis lineal
vGene <- voom(dge, mod, plot = TRUE)

# Ajuste del modelo lineal y cálculo de estadísticas empíricas de Bayes
eb_results <- eBayes(lmFit(vGene))


# Extraer la tabla de genes diferencialmente expresados.
de_results <- topTable(
  eb_results,
  coef = 2, # Se refiere al coeficiente del segundo término en el modelo
  number = nrow(rse_gene_SRP075398),
  sort.by = "none"
)

# Dimensiones y vista preliminar de los resultados
dim(de_results)

head(de_results)

# Genes diferencialmente expresados con FDR < 5%
table(de_results$adj.P.Val < 0.05)

# Visualizar los resultados estadísticos

plotMA(eb_results, coef = 2)

volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

# Información de los 3 genes más significativos
de_results[de_results$gene_name %in% c("AMIGO2", "AFAP1L2", "PHLDA1"), ]

# Revisar los top 50 genes diferencialmente expresados

# Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

# Crear una tabla con información de las muestras y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP075398)[, c("sra_attribute.cell_line", "sra_attribute.transfection")])
colnames(df) <- c("Cell_line", "Transfection")

# Hacer un heatmap
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

## Guardemos los IDs de nuestros 50 genes
nombres_originales <- rownames(exprs_heatmap)

## Con match() podemos encontrar cual es cual
rownames(exprs_heatmap) <- rowRanges(rse_gene_SRP075398)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(rse_gene_SRP075398)$gene_id)
]

## Guardar la imagen en un PDF largo para poder ver los nombres de los genes
pdf("pheatmap_con_nombres.pdf", height = 16, useDingbats = FALSE)
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  fontsize_row = 6,
)

dev.off()

# MDS (multidimensional scaling)

## Convertir los grupos de Cell_line a colores
col.group <- df$Cell_line
levels(col.group) <- brewer.pal(nlevels(col.group), "Set2")
col.group <- as.character(col.group)

## MDS por grupos de Cell_line
plotMDS(vGene$E, labels = df$Cell_line, col = col.group)

## Convertir Transfection a colores
col.group <- df$Transfection
df$Transfection <- as.factor(df$Transfection) # Asegúrate de que sea un factor
colors <- brewer.pal(nlevels(df$Transfection), "Set2") # Generar paleta de colores
levels(col.group) <- colors # Asignar colores a los niveles
col.group <- as.character(col.group) # Convertir a vector de caracteres

## MDS por grupos de Transfection
plotMDS(vGene$E, labels = df$Transfection, col = col.group,
        main = "MDS Plot by Transfection", pch = 16)

## Agregar leyenda
legend("topright", legend = levels(df$Transfection), fill = unique(col.group),
       title = "Transfection Groups")

