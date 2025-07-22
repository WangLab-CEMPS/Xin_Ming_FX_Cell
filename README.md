# FX-Cell: a method for single-cell RNA sequencing on difficult-to-digest and cryopreserved plant samples

The FX-Cell method and its derivative technologies, FXcryo-Cell and cryoFX-Cell, were developed to address the challenges of processing plant samples that are difficult to digest and cryopreserve for single-cell RNA sequencing (scRNA-seq). By overcoming these obstacles, these methods have significantly expanded the application of scRNA-seq in plants, yielding high-quality single-cell data from a variety of tissues and organs and thereby promoting the widespread application of plant single-cell genomics research.

> ***This repository contains the code for the downstream analysis of FX-Cell data.***

## Document Structure

- `fig_scripts`: The code used to generate the majority of the main and supplementary figures.

- `src`: Specific analysis code for each of the species studied.
    - `athleaf`: Arabidopsis thaliana leaf
        - `snRNA`: Single-nucleus RNA-seq
    - `athroot`: Arabidopsis thaliana root
    - `maize`: Maize crown root
    - `osroot`: Oryza sativa root
        - `ggm`: FX-Cell
        - `dgm`: cryoFX-Cell
        - `gdm`: FXcryo-Cell
    - `tillering_node_rhizome`: Rice tiller nodes and wild rice rhizome nodes
    - `origin_method`: snRNA-seq analysis of Oryza sativa root