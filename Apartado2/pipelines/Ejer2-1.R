################################################################################
################################ CARGAR PAQUETES ###############################
################################################################################

library("RColorBrewer")
library("genefilter")
library("gplots")
library("ggplot2")
library("ggrepel")
library("rlang")
library("dplyr")
library("readr")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("DESeq2")
library("pheatmap")


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

# 1.7 SUBSET DE LOS RESULTADOS A 24h
dds_24 <- dds[, dds$time == "24h"]
colData(dds_24)

################################################################################
#################### 2. ANÁLISIS EXPLORATORIO DE LOS DATOS #####################
################################################################################

# ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)

vsd <- vst(dds_24, blind=TRUE)
print(plotPCA(vsd, intgroup = c("patient", "agent")))

#### Conclusión: eliminamos al paciente 4

# ANÁLISIS DE DISTANCIAS

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$patient, vsd$agent,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clistering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)

################################################################################
########### 3. PREPARAR EL OBJETO PARA EL ANÁLISIS Y CORRER PIPELINE ###########
################################################################################

# 3.1 ELIMINAR AL PACIENTE 4 (OUTLIER) 
dds_24_filt <- dds_24[, dds_24$patient != "4"] 
dds_24_filt$patient <- factor(dds_24_filt$patient) #factorizar de nuevo si es necesario
colData(dds_24_filt)

vsd_filt <- vst(dds_24_filt, blind=TRUE)
print(plotPCA(vsd_filt, intgroup = c("patient", "agent")))

sampleDists <- dist(t(assay(vsd_filt)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_filt$patient, vsd_filt$agent,
                                    sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clistering_distance_rows=sampleDists, 
         clustering_distance_cols=sampleDists, 
         col=colors)


# 3.2 CORRER EL PIPELINE E INSPECCIONAR LOS RESULTADOS A LAS 24h
dds_24_filt <- DESeq(dds_24_filt, test = "Wald")
res_24_filt <- results(dds_24_filt)
res_24_filt
summary(res_24_filt)
plotDispEsts(dds_24_filt)

# (OPCIONAL) añadimos los nombres de los genes asociados a cada símbolo
res_24_filt$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_24_filt), 
                             keytype = "ENSEMBL", column = "SYMBOL")
res_24_filt

#### 3.3 GENES DIFERENCIALMENTE EXPRESADOS TRAS 24H DE TRATAMIENTO CON DPN ####

#### obtenemos los DEGs para el tratamiento con DPN
res_DPN_24 <- results(dds_24_filt, contrast = c("agent", "DPN", "Control"), 
                      alpha=0.05)
res_DPN_24$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_DPN_24), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res_DPN_24
summary(res_DPN_24) #número de genes up y downregualados

#### PREPARAMOS LOS DATOS PARA EL VOLCANO PLOT
DPN_DEGs_24 <- results(object = dds_24_filt,
                       contrast = c("agent", "DPN", "Control"),
                       alpha = 0.05,
                       pAdjustMethod = "BH",
                       tidy = TRUE
)

DPN_DEGs_24$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_24_filt), 
                             keytype = "ENSEMBL", column = "SYMBOL")

DPN_DEGs_24


## definimos los DEGs como aquellos significativos, siendo "NO" los no significativos
DPN_DEGs_24$diffexpressed <- "NO"
DPN_DEGs_24$diffexpressed[DPN_DEGs_24$log2FoldChange > 0  & 
                            DPN_DEGs_24$padj < 0.05] <- "UP"
DPN_DEGs_24$diffexpressed[DPN_DEGs_24$log2FoldChange < -0 & 
                            DPN_DEGs_24$padj < 0.05] <- "DOWN"

DPN_DEGs_24[DPN_DEGs_24$diffexpressed != "NO", ] 

#### definimos la funcion delabel para asignar el genesymbol solo a los signif.
DPN_DEGs_24$delabel <- NA
DPN_DEGs_24$delabel[DPN_DEGs_24$diffexpressed !="NO"] <- 
  DPN_DEGs_24$symbol[DPN_DEGs_24$diffexpressed !="NO"]

#### visualizamos el volcano plot: en rojo los urpregulados y en azul los downregulados
ggplot(data=DPN_DEGs_24, aes(x=log2FoldChange, y=-log10(pvalue), 
                             col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=c("blue", "black", "red"))
  

#### MA PLOT
plotMA(res_DPN_24)

#### idem estableciendo lfc=1
DPN_DEGs_lfc1 <- results(object = dds_24_filt,
                         contrast = c("agent", "DPN", "Control"),
                         lfcThreshold = 1,
                         alpha = 0.05,
                         pAdjustMethod = "BH",
)


DPN_DEGs_lfc1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(dds_24_filt), 
                               keytype = "ENSEMBL", column = "SYMBOL")

summary(DPN_DEGs_lfc1)

#### 3.3 GENES DIFERENCIALMENTE EXPRESADOS TRAS 24H DE TRATAMIENTO CON OHT ####

#### obtenemos los DEGs para el tratamiento con OHT
res_OHT_24 <- results(dds_24_filt, contrast = c("agent", "OHT", "Control"), 
                      alpha=0.05)
res_OHT_24$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_OHT_24), 
                            keytype = "ENSEMBL", column = "SYMBOL")
res_OHT_24
summary(res_OHT_24) #número de genes up y downregualados

#### PREPARAMOS LOS DATOS PARA EL VOLCANO PLOT
OHT_DEGs_24 <- results(object = dds_24_filt,
                       contrast = c("agent", "OHT", "Control"),
                       alpha = 0.05,
                       pAdjustMethod = "BH",
                       tidy = TRUE
)

OHT_DEGs_24$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_24_filt), 
                             keytype = "ENSEMBL", column = "SYMBOL")

OHT_DEGs_24


## definimos los DEGs como aquellos significativos, siendo "NO" los no significativos
OHT_DEGs_24$diffexpressed <- "NO"
OHT_DEGs_24$diffexpressed[OHT_DEGs_24$log2FoldChange > 0  & 
                            OHT_DEGs_24$padj < 0.05] <- "UP"
OHT_DEGs_24$diffexpressed[OHT_DEGs_24$log2FoldChange < -0 & 
                            OHT_DEGs_24$padj < 0.05] <- "DOWN"

OHT_DEGs_24[OHT_DEGs_24$diffexpressed != "NO", ] 

#### definimos la funcion delabel para asignar el genesymbol solo a los signif.
OHT_DEGs_24$delabel <- NA
OHT_DEGs_24$delabel[OHT_DEGs_24$diffexpressed !="NO"] <- 
  OHT_DEGs_24$symbol[OHT_DEGs_24$diffexpressed !="NO"]

#### visualizamos el volcano plot: en rojo los urpregulados y en azul los downregulados
ggplot(data=OHT_DEGs_24, aes(x=log2FoldChange, y=-log10(pvalue), 
                             col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=c("black", "blue", "red"))

#### MA PLOT
plotMA(res_OHT_24)

#### idem estableciendo lfc=1
OHT_DEGs_lfc1 <- results(object = dds_24_filt,
                         contrast = c("agent", "OHT", "Control"),
                         lfcThreshold = 1,
                         alpha = 0.05,
                         pAdjustMethod = "BH"
)


OHT_DEGs_lfc1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(dds_24_filt), 
                               keytype = "ENSEMBL", column = "SYMBOL")

summary(OHT_DEGs_lfc1)
