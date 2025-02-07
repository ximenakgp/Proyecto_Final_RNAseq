# CÓDIGO PARA DETERMINAR EL MODELO Y VISUALIZAR LA MATRIZ

library(ExploreModelMatrix)
library(cowplot)

# Construcción de la matriz de diseño para el modelo lineal.
mod <- model.matrix(
  ~ sra_attribute.cell_line + sra_attribute.transfection + assigned_gene_prop,
  data = colData(rse_gene_SRP075398)
)

# Cada columna representa un coeficiente del modelo
colnames(mod)

# Visualizar la matriz de diseño completa
mod

## Crear las visualizaciones
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = colData(rse_gene_SRP075398), # Metadatos de las muestras
  designFormula = ~ sra_attribute.cell_line + sra_attribute.transfection,
  textSizeFitted = 2
)

cowplot::plot_grid(plotlist = vd$plotlist)


