# Tables and Figure for CN multiome
# LZ 28.10.2024


library(readr)
library(readxl)
library(data.table)


setwd("/path/to/")


#### DE Supplementary Table S2

#read results
astro_DE <- read_csv("results/DE/FindMarkers/all/clusterAstrocyte.csv")
micro_DE <- read_csv("results/DE/FindMarkers/all/clusterMicroglia.csv")
oligo_DE <- read_csv("results/DE/FindMarkers/all/clusterOligodendrocyte.csv")
OPC_DE <- read_csv("results/DE/FindMarkers/all/clusterOPC.csv")
D1_DE <- read_csv("results/DE/FindMarkers/all/clusterGABAergic - D1 MSN.csv")
D2_DE <- read_csv("results/DE/FindMarkers/all/clusterGABAergic - D2 MSN.csv")

colnames(astro_DE) <- paste0(colnames(astro_DE),"_astro")
colnames(micro_DE) <- paste0(colnames(micro_DE),"_micro")
colnames(oligo_DE) <- paste0(colnames(oligo_DE),"_oligo")
colnames(OPC_DE) <- paste0(colnames(OPC_DE),"_OPC")
colnames(D1_DE) <- paste0(colnames(D1_DE),"_D1")
colnames(D2_DE) <- paste0(colnames(D2_DE),"_D2")

astro_DE <- astro_DE[order(astro_DE$...1_astro),]
micro_DE <- micro_DE[order(micro_DE$...1_micro),]
oligo_DE <- oligo_DE[order(oligo_DE$...1_oligo),]
OPC_DE <- OPC_DE[order(OPC_DE$...1_OPC),]
D1_DE <- D1_DE[order(D1_DE$...1_D1),]
D2_DE <- D2_DE[order(D2_DE$...1_D2),]

#combine dataframes

DEall <- cbind(astro_DE, D1_DE, D2_DE, micro_DE, oligo_DE, OPC_DE)
DEall_sig <- DEall[DEall$p_val_adj_astro<0.05|DEall$p_val_adj_micro<0.05|DEall$p_val_adj_oligo<0.05|
                     DEall$p_val_adj_OPC<0.05|DEall$p_val_adj_D1<0.05|DEall$p_val_adj_D2<0.05,]

DEall_sig <- DEall_sig[order(DEall_sig$p_val_adj_astro),]

write_csv(DEall_sig, "results/DE/FindMarkers/DE_all_cell_types.csv")



#### DA Supplementary Table S3
#read results
astro_DA <- read_csv("results/DAc/FindMarkers/clusterAstrocyte.csv")
micro_DA <- read_csv("results/DAc/FindMarkers/clusterMicroglia.csv")
oligo_DA <- read_csv("results/DAc/FindMarkers/clusterOligodendrocyte.csv")
OPC_DA <- read_csv("results/DAc/FindMarkers/clusterOPC.csv")
D1_DA <- read_csv("results/DAc/FindMarkers/clusterGABAergic - D1 MSN.csv")
D2_DA <- read_csv("results/DAc/FindMarkers/clusterGABAergic - D2 MSN.csv")

colnames(astro_DA) <- paste0(colnames(astro_DA),"_astro")
colnames(micro_DA) <- paste0(colnames(micro_DA),"_micro")
colnames(oligo_DA) <- paste0(colnames(oligo_DA),"_oligo")
colnames(OPC_DA) <- paste0(colnames(OPC_DA),"_OPC")
colnames(D1_DA) <- paste0(colnames(D1_DA),"_D1")
colnames(D2_DA) <- paste0(colnames(D2_DA),"_D2")

astro_DA <- astro_DA[order(astro_DA$...1_astro),]
micro_DA <- micro_DA[order(micro_DA$...1_micro),]
oligo_DA <- oligo_DA[order(oligo_DA$...1_oligo),]
OPC_DA <- OPC_DA[order(OPC_DA$...1_OPC),]
D1_DA <- D1_DA[order(D1_DA$...1_D1),]
D2_DA <- D2_DA[order(D2_DA$...1_D2),]

#combine dataframes

DAall <- cbind(astro_DA, D1_DA, D2_DA, micro_DA, oligo_DA, OPC_DA)
DAall_sig <- DAall[DAall$p_val_adj_astro<0.05|DAall$p_val_adj_micro<0.05|DAall$p_val_adj_oligo<0.05|
                     DAall$p_val_adj_OPC<0.05|DAall$p_val_adj_D1<0.05|DAall$p_val_adj_D2<0.05,]

DAall_sig <- DAall_sig[order(DAall_sig$p_val_adj_astro),]

write_csv(DAall_sig, "results/DAc/FindMarkers/DA_all_cell_types.csv")


#### Motifs Supplementary Table S4
#read results
astro_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_Astrocyte.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
micro_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_Microglia.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
oligo_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_Oligodendrocyte.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
OPC_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_OPC.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
D1_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_GABAergic - D1 MSN.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
D2_motif <- read_delim("results/chromVar/Motif_Enrichment_cluster_GABAergic - D2 MSN.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)

colnames(astro_motif) <- paste0(colnames(astro_motif),"_astro")
colnames(micro_motif) <- paste0(colnames(micro_motif),"_micro")
colnames(oligo_motif) <- paste0(colnames(oligo_motif),"_oligo")
colnames(OPC_motif) <- paste0(colnames(OPC_motif),"_OPC")
colnames(D1_motif) <- paste0(colnames(D1_motif),"_D1")
colnames(D2_motif) <- paste0(colnames(D2_motif),"_D2")

astro_motif <- astro_motif[order(astro_motif$motif_astro),]
micro_motif <- micro_motif[order(micro_motif$motif_micro),]
oligo_motif <- oligo_motif[order(oligo_motif$motif_oligo),]
OPC_motif <- OPC_motif[order(OPC_motif$motif_OPC),]
D1_motif <- D1_motif[order(D1_motif$motif_D1),]
D2_motif <- D2_motif[order(D2_motif$motif_D2),]

#combine dataframes

motifall <- cbind(astro_motif, D1_motif, D2_motif, micro_motif, oligo_motif, OPC_motif)
motifall_sig <- motifall[motifall$p_val_adj_astro<0.05|motifall$p_val_adj_micro<0.05|motifall$p_val_adj_oligo<0.05|
                     motifall$p_val_adj_OPC<0.05|motifall$p_val_adj_D1<0.05|motifall$p_val_adj_D2<0.05,]

motifall_sig <- motifall_sig[order(motifall_sig$p_val_adj_astro),]

write_csv(motifall_sig, "results/chromVar/motifs_all_cell_types.csv")
