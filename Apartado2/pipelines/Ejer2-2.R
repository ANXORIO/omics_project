################################################################################
################################ CARGAR PAQUETES ###############################
################################################################################
rm(list = ls()) # R version 4.0.3

library("stringr") # stringr_1.4.0
library("apeglm")
library("RColorBrewer")
library("genefilter")
library("gplots")
library("DESeq2")
library("ggplot2")
library("dplyr")
library("readr")
library("AnnotationDbi")
library("org.Hs.eg.db")

################################################################################
############################### DEFINIR FUNCIONES ##############################
################################################################################

output.gmt <- function(l, filename) {
  if(file.exists(filename)) {
    warning(paste("Removing previous", filename, "file."))
    file.remove(filename)
  }
  invisible(lapply(1:length(l), function(i) {
    cat(c(names(l)[i], "na", l[[i]], "\n"), sep = "\t", file = filename,
        append = TRUE)
  }))
}

################################################################################
#### 1.GENERAMOS EL DATASET DE TRABAJO A PARTIR DE LAS RAWCOUNTS Y METADATA ####
################################################################################

# 1.1 FIJAR EL DIRECTORIO DE TRABAJO DONDE SE ENCUENTRAN LOS ARCHIVOS

WD <- setwd("/home/vant/transcriptomics-project/Apartado2/input")
WD

# 1.2 CARGAR LAS TABLAS RAWCOUNTS Y METADATA

metadata <- read.table(file="metadata.tsv")
rawcounts <- read.table(file="rawcounts.tsv")

# 1.3 FIJAR LOS PARÁMETROS DEL METADATA COMO FACTORES

metadata$patient <- factor(metadata$patient)
metadata$agent <- factor(metadata$agent)
metadata$time <- factor(metadata$time)

# 1.4 COMPROBAR COINCIDENCIA ENTRE COLUMNAS DE COUNTS Y FILAS DE METADATA

all(colnames(rawcounts) == rownames(metadata))
all(colnames(rawcounts) %in% rownames(metadata))

#### continuar con el analisis si el resultado es TRUE ####

# 1.5 GENERAR EL DATASET A PARTIR DE LA MATRIZ DE DATOS

dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = 
                                metadata, design = ~ patient + agent)

# 1.6 ELIMINAR TODAS LAS FILAS QUE NO TENGAN AL MENOS 10 LECTURAS

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# 1.7 SUBSET DE LOS RESULTADOS A 24h Y 48h
dds_24 <- dds[, dds$time == "24h"]
colData(dds_24)

dds_48 <- dds[, dds$time == "48h"]
colData(dds_48)

################################################################################
#################### 2. ANÁLISIS EXPLORATORIO DE LOS DATOS #####################
################################################################################

# ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)

vsd_24 <- vst(dds_24, blind=TRUE)
print(plotPCA(vsd_24, intgroup = c("patient", "agent")))

vsd_48 <- vst(dds_48, blind=TRUE)
print(plotPCA(vsd_48, intgroup = c("patient", "agent")))


#### Conclusión: eliminamos al paciente 4 en la condición 24h

# ANÁLISIS DE DISTANCIAS EN TRATAMIENTO 24H

sampleDists <- dist(t(assay(vsd_24)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_24$patient, vsd$agent,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clistering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)

# ANÁLISIS DE DISTANCIAS EN TRATAMIENTO 48H

sampleDists <- dist(t(assay(vsd_48)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_48$patient, vsd$agent,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clistering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)

################################################################################
#################### 3. PREPARAR EL OBJETO PARA EL ANÁLISIS  ###################
################################################################################

# 3.1 ELIMINAR AL PACIENTE 4 (OUTLIER) EN LA CONDICIÓN 24h
dds_24_filt <- dds_24[, dds_24$patient != "4"] 
dds_24_filt$patient <- factor(dds_24_filt$patient) #factorizar de nuevo si es necesario
colData(dds_24_filt)

vsd_24_filt <- vst(dds_24_filt, blind=TRUE)
print(plotPCA(vsd_24_filt, intgroup = c("patient", "agent")))

sampleDists <- dist(t(assay(vsd_24_filt)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_24_filt$patient, vsd_24_filt$agent,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clistering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)


#### 3.2.1 CORRER EL PIPELINE E INSPECCIONAR LOS RESULTADOS A LAS 24h
dds_24_filt <- DESeq(dds_24_filt, test = "Wald")
res_24_filt <- results(dds_24_filt)
res_24_filt
summary(res_24_filt)
plotDispEsts(dds_24_filt)

#### 3.2.1 CORRER EL PIPELINE E INSPECCIONAR LOS RESULTADOS A LAS 48h
dds_48 <- DESeq(dds_48, test = "Wald")
res_48 <- results(dds_48)
res_48
summary(res_48)
plotDispEsts(dds_48)

# (OPCIONAL) añadimos los nombres de los genes asociados a cada símbolo
res_24_filt$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_24_filt), 
                             keytype = "ENSEMBL", column = "SYMBOL")
res_24_filt

res_48$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_48), 
                             keytype = "ENSEMBL", column = "SYMBOL")
res_48


################################################################################
########### 4. GENES DIFERENCIALMENTE EXPRESADOS TRATAMIENTO CON DPN ###########
################################################################################

#### TOP 100 GENES UP/DOWNREGULADOS TRAS 24 HORAS DE TRATAMIENTO

res_DPN_24 <- results(dds_24_filt, contrast = c("agent", "DPN", "Control"), 
                      alpha=0.05, pAdjustMethod = "BH")
res_DPN_24$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_DPN_24), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res_DPN_24
summary(res_DPN_24) #número de genes up y downregulados

resultsNames(dds_24_filt)
res_ape_24 <- lfcShrink(dds_24_filt, coef = "agent_DPN_vs_Control", type = "apeglm",
                        res = res_DPN_24)
summary(res_ape_24) # Same number of up/down genes padj < 0.05

rnk_24 <- data.frame(Feature = rownames(res_ape_24), LFC = res_ape_24$log2FoldChange)
head(rnk_24)
rnk_24$Feature <- str_remove(rnk_24$Feature, "\\..*$")
head(rnk_24)

write.table(rnk_24, file = "DEGs_24.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)


#### FIRMAS TOP 100 GENES UP/DOWNREGULADOS TRAS 48 HORAS DE TRATAMIENTO
res_DPN_48 <- results(dds_48, contrast = c("agent", "DPN", "Control"), 
                      alpha=0.05, pAdjustMethod = "BH")
res_DPN_48$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_DPN_48), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res_DPN_48
summary(res_DPN_48) #número de genes up y downregualados

resultsNames(dds_48)
res_ape_48 <- lfcShrink(dds_48, coef = "agent_DPN_vs_Control", type = "apeglm",
                        res = res_DPN_48)
summary(res_ape_48) # Same number of up/down genes padj < 0.05

rnk_48 <- data.frame(Feature = rownames(res_ape_48), LFC = res_ape_48$log2FoldChange)
head(rnk_48)
rnk_48$Feature <- str_remove(rnk_48$Feature, "\\..*$")
head(rnk_48)

write.table(rnk_48, file = "DEGs_48.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)

up100_48 <- rnk_48 %>% arrange(desc(LFC))
head(up100_48)
up100_48 <- up100_48$Feature[1:100]

down100_48 <- rnk_48 %>% arrange(LFC)
head(down100_48)
down100_48 <- down100_48$Feature[1:100]

geneset_48 <- list(up_48 = up100_48, down_48 = down100_48)
output.gmt(geneset_48, filename = "DEGs_48.gmt")
