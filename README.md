# A Hitchhiker's guide to high-dimensional tissue imaging with multiplexed ion beam imaging

## Abstract
Advancements in multiplexed tissue imaging technologies are vital in shaping our understanding of tissue microenvironmental influences in disease contexts. These technologies now allow us to relate the phenotype of individual cells to their higher-order roles in tissue organization and function. Multiplexed Ion Beam Imaging (MIBI) is one of such technologies, which uses metal isotope-labeled antibodies and secondary ion mass spectrometry (SIMS) to image more than 40 protein markers simultaneously within a single tissue section. Here, we describe an optimized MIBI workflow for high-plex analysis of Formalin-Fixed Paraffin-Embedded (FFPE) tissues following antigen retrieval, metal isotope-conjugated antibody staining, imaging using the MIBI instrument, and subsequent data processing and analysis. While this workflow is focused on imaging human FFPE samples using the MIBI, this workflow can be easily extended to model systems, biological questions, and multiplexed imaging modalities.

## Overview of code
The following is an overview of the code used in this repository. These are stored under "Scripts".
| File Name | Description |
| :---------- | :---------- |
| 1a. prepare_tif.ijm | Example of an ImageJ macro that helps to normalize bit size and perform signal capping on raw images
| 1b. prepare_surface.ijm | Example of an ImageJ macro to generate a cell surface marker
| 2a. Mesmer_custom | Jupyter notebook for cell segmentation
| 2b. deepcell_check.ijm | ImageJ macro to help visualize segmentation by overlaying with nucleus and surface markers
| 3. extractCellData.m | Single-cell feature extraction, based on cell segmentation
| 3. writeFCS.m | Contains helper functions for the extractCellData.m script
| 4. cellClustering.rmd | Signal capping and normalization, cell clustering, visualization of clusters by UMAP, and visualization of cluster features MEM heatmap
| 4. MIBI_functions.R | Contains helper function for the cellClustering.rmd script
| 5a. cellColoring.m | Example code for coloring cells based on segmentation mask
| 5b. cellColoringCluster.m | Produce cell cluster maps by coloring each cluster; this helps to visually verify annotated cell types
| 5c. cellColoringClusterAnnotate.m | Produce phenotype map by coloring annotated cell types

## Sample data
The following is an overview of example data generated as part of the figure for the manuscript. These are stored under "Scripts".
| File Name | Description |
| :---------- | :---------- |
| Annotate_output | Generated phenotype map, as well as cell-type specific phenotype maps
| Cell_segmentation | Cell nuclei and surface images, as well as resulting segmentation outlines and masks
| Cluster_output | Cell cluster maps to assist in visual verification of annotated cell types
| FCS_output | Extracted single-cell features for cell clustering including cell-size scaled and cell-size unscaled versions
| TIFsNoAgg | Raw MIBI images (with background removed, denoised, and aggregates removed)