# Ajustar las variables al formato correcto para usarlas en el modelo estadístico

## Convertir las variables categóricas de character a factor

colData(rse_gene_SRP068565)$drug_treatment <- factor(colData(rse_gene_SRP068565)$sra_attribute.drug_treatment)

colData(rse_gene_SRP068565)$source_name <- factor(colData(rse_gene_SRP068565)$sra_attribute.source_name)

## Resumen de las variables

summary(as.data.frame(colData(rse_gene_SRP068565)[
  ,
  grepl("^sra_attribute.[drug_treatment|source_name]", colnames(colData(rse_gene_SRP068565)))
]))

# sra_attribute.drug_treatment sra_attribute.source_name
# Length:20                    Length:20
# Class :character             Class :character
# Mode  :character             Mode  :character

