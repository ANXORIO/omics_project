library("DESeq2")

# DEFINIR EL DIRECTORIO DE TRABAJO
wd <- setwd("/home/vant/transcriptomics-project/Apartado1/input/htseq")

# ALMACENAR LOS NOMBRES DE LAS MUESTRAS Y ASIGNAR CONDICIONES
sampleFiles <- grep("htseq", list.files(wd), value=TRUE)
sampleCondition <- c("control", "control", "treated", "treated")
sampleTable <- data.frame(sampleName = sampleFiles, fileName = 
                            sampleFiles, condition = sampleCondition)
# ESTABLECER LAS CONDICIONES COMO FACTORES
sampleTable$condition <- factor(sampleTable$condition)

# GENERAR LA MATRIZ DE DATOS
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = wd,
                                  design= ~ condition)
dds

# ELIMINAR VALORES CON MENOS DE 10 COUNTS
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# EXPORTAR LOS DATOS DE COUNTS CRUDOS
rawcounts <- counts(dds)
write.csv(rawcounts, file="/home/vant/transcriptomics-project/Apartado1/output/rawcounts.csv")

# EXTRAER LOS DATOS DE COUNTS NORMALIZADOS
dds <- DESeq(dds)
normcounts <- counts(dds, normalized=TRUE)
write.csv(normcounts, file="/home/vant/transcriptomics-project/Apartado1/output/normcounts.csv")
