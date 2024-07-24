library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)
library(ggplot2)
library(knitr)
set.seed(42)

for (i in c(paste0("sample",1:16))){

counts <- Read10X_h5(paste0("/path/to//cellranger_out/",i,"/outs/filtered_feature_bc_matrix.h5"))

# Annotations 
ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 105 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"

# Create Seurat object containing RNA and ATAC data
sample <- CreateSeuratObject(counts = counts$`Gene Expression`,
                               assay = "RNA",project = i)
sample[['ATAC']] <- CreateChromatinAssay(counts = counts$`Peaks`,
                                           annotation = annotations,
                                           fragments = paste0("/zi-flstorage/group_genepi/data/single-cell_data/2023_03_24_CUD_Multiome/cellranger_out/",i,"/outs/atac_fragments.tsv.gz"),sep = c(":", "-"),genome = 'hg38')

sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
sample <- NucleosomeSignal(sample, assay = "ATAC")
sample <- TSSEnrichment(sample, assay = "ATAC")

sample <- subset(sample,
                   subset = nFeature_RNA > 200 &
                     nFeature_RNA < 6500 &
                     percent.mt < 10 &
                     nFeature_ATAC > 1000 &
                     nFeature_ATAC < 25000 &
                     TSS.enrichment > 2 &
                     TSS.enrichment < 8 &
                     nucleosome_signal < 2 &
                     nucleosome_signal > 0.5)

saveRDS(sample,file = paste0("/path/to/data/postQC/",i,".rds"))
}
