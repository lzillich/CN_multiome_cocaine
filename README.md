## Cell type-specific Multi-Omics Analysis of Cocaine Use Disorder in the Human Caudate Nucleus (Zillich et al., 2025)

* 1_initial_QC - QC scripts per sample 
* 2_merge_call_peaks - seurat files are merged and a joint peak calling is performed (MACS2)
* 3a_RNA - Preprocessing of RNA slot
* 3b_ATAC_RNA_integration - Preprocessing of ATAC slot, integration with RNA slot using harmony
* 4_Markers_Integrated_Figure1_S5 - Cluster annotation, code to produce Figure 1 and Supplementary Figure 5
* 5a_MOFA_prep - create MOFA object, run MOFA
* 5b_MOFA_downstream_FigureS4 - downstream analysis of MOFA, code to produce Supplementary Figure 4A-4E
* 6_DE_DA_analysis_Figure2_FigureS2 - differential expression and accessibility analysis, code to produce Figures 2 and Supplementary Figure 2
* 7_motifs_TF_Figure3 - motif analyses, code to produce Figure 3
* 8_Pando_MSN_Figure4_S3 - pando in D1/D2 medium spiny neurons, code to produce Figure 4A-4C and Supplementary Figure 3
* 9a_scDRS_FigureS4 - scDRS analysis
* 9b_visualize_scDRS_FigureS4 - visualisation of scDRS results, code to produce Supplementary Figure 4F 
* 10_RRHO_sex_specific.Rmd - sex-specific DE and RRHO analysis, code to produce Figure 4E
* 11_qPCR.R - analysis of qPCR results, code to produce Figure 4D
* 12_tables_manuscript.R - code to produce Supplementary Tables
