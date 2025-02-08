## Proyecto de Análisis de Expresión Diferencial

Este repositorio contiene los archivos necesarios para el análisis de expresión diferencial a partir de datos de secuenciación de RNA (RNA-seq) obtenidos con recount3 del estudio SRP075398. 
Incluye el código en R, el documento en formato R Markdown, las figuras generadas, y los resultados en diversos formatos.

## Estructura del repositorio
1. **Directorios**

  - Proyecto_RNAseq_files/figure-latex: Contiene las figuras generadas en formato LaTeX durante la ejecución del archivo R Markdown.
  
  - R/: Contiene scripts adicionales en R utilizados en el análisis.

2. **Archivos principales**

  - .gitignore: Lista de archivos y carpetas que Git debe ignorar en el repositorio.
  
  - Proyecto_Final_RNAseq.Rproj: Archivo de proyecto de RStudio, que facilita la organización y gestión del entorno de trabajo.
  
  - Proyecto_RNAseq.Rmd: Archivo principal en formato R Markdown que documenta el análisis de expresión diferencial, incluyendo código R, resultados y visualizaciones.
  
  - Proyecto_RNAseq.html: Versión renderizada en formato HTML del documento R Markdown, mostrando los resultados y figuras generadas.
  
  - Proyecto_RNAseq.pdf: Versión en formato PDF del análisis generado a partir del archivo R Markdown.
  
  - Proyecto_RNAseq.tex: Archivo en formato LaTeX generado a partir del R Markdown, utilizado para la conversión a PDF.
  
  - export_references.bib: Archivo en formato BibTeX que contiene referencias bibliográficas utilizadas en el análisis.
  
  - pheatmap_con_nombres.pdf: Gráfico de mapa de calor (heatmap) con nombres de genes o condiciones, generado como parte del análisis de expresión diferencial.

## Uso

Para reproducir el análisis, abre el archivo Proyecto_RNAseq.Rmd en RStudio y ejecúta el código paso a paso o renderiza el documento completo con:

``` R
rmarkdown::render("Proyecto_RNAseq.Rmd")
```
Este comando generará los resultados en los formatos HTML, PDF o LaTeX según la configuración del archivo.

## Sobre el Proyecto

Este proyecto fue desarrollado como parte del módulo de Introducción a R, RStudio y Secuenciación de ARN (RNA-seq), impartido en enero de 2025 para el 
programa LCG-UNAM en el CCG-UNAM (del 28 al 31 de enero de 2025). 
El curso fue dirigido por el Dr. Leonardo Collado-Torres, cuyo material puede encontrarse en el siguiente repositorio: <https://lcolladotor.github.io/rnaseq_LCG-UNAM_2025/>

## Autor

Este repositorio fue creado y mantenido por ximenakgp para el análisis de expresión diferencial en RNA-seq.
Contáctenos: [ximenagp@lcg.unam.mx]
