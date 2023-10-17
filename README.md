## Description
Multiome data (RNA-seq + ATAC-seq) is composed of two cachectic replicates and two controls (SHAM). For each replicate, we create a Seurat object with RNA and ATAC. 
- Quality control (1_multiome_quality_control.R);
- RNA and ATAC integration and doublets removal (2_multiome_RNA_ATAC_integration.R);
- clustering and cell-type annotation (3_multiome_celltype_annotation.R);
- differential expression analysis and differential abundance (4_multiome_DEA_DA.R);

GEX data is composed of one cachectic and sham replicate. For each replicate, we create a Seurat object.
- Quality control (5_rna_quality_control.R );
- RNA integration (6_rna_integration.R);
- clustering and cell-type annotation (7_rna_celltype_annotation.R);

Multiome + RNA 
- multiome and rna data integration (8_multiome_rna_integration.R);
- clustering and cell-type annotation (9_multiome_rna_celltype_annotation.R);
- differential expression analysis and differential abundance (10_multiome_rna_DEA_DA.R).

