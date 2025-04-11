## Description
Multiome data (RNA-seq + ATAC-seq) is composed of two cachectic replicates and two controls (SHAM). For each replicate, we create a Seurat object with RNA and ATAC. 
- Quality control (1_quality_control.R);
- RNA and ATAC integration and doublets removal (2_rna_atac_integration.R);
- clustering and cell-type annotation (3_celltype_annotation.R);
- differential expression analysis and differential abundance (4_DEA_DAA.R);


