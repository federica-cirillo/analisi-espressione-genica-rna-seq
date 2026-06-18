# RNA Sequencing Analysis of GSE267074

A complete bulk RNA sequencing analysis of the public GEO dataset GSE267074. The data come from an experiment on immunocompromised mice carrying xenografts of human B cell non Hodgkin lymphoma, an aggressive cancer of the lymphatic system. The mice were treated with CAR T cells directed against CD19, given either together with tazemetostat (a selective inhibitor of the epigenetic regulator EZH2) or with a vehicle control. The dataset contains 16,282 genes measured across 18 samples, corresponding to six mice (three per treatment group) with three technical replicates each. The goal is to uncover the transcriptional changes and the biological pathways driven by the combined CAR T and tazemetostat treatment, in order to clarify the molecular basis of a possible synergy between epigenetic inhibition and CAR T immunotherapy.

## Workflow

The analysis starts with data preprocessing. Raw counts and sample metadata are downloaded directly from GEO and assembled into a SummarizedExperiment object, lowly expressed genes are filtered out, and data quality is inspected through boxplots, RLE plots and PCA. Several normalization strategies are compared (upper quantile and full quantile from EDASeq, TMM and RLE from edgeR), and one outlier sample (GSM8259456) is detected and removed, with its impact assessed by repeating the analysis with and without it.

Differential expression is then performed with edgeR, fitting a generalized linear model and testing every gene with a likelihood ratio test on the treatment contrast. The significant genes are explored through heatmaps, p value histograms and MA plots.

Finally, the differentially expressed genes are studied through functional enrichment analysis. Gene Ontology terms across all three domains (biological process, molecular function and cellular component) and KEGG pathways are tested with clusterProfiler and visualized through dotplots, enrichment maps, gene concept networks and detailed pathway graphs, surfacing immune related signaling such as the B cell receptor and Toll like receptor pathways.

## Repository structure

The report is provided both as an R Markdown source (`rnaseq_analysis_GSE267074.Rmd`) and as its rendered HTML version (`rnaseq_analysis_GSE267074.html`), which opens in any browser and shows the full analysis with code, figures and comments. The custom plotting and helper functions live separately in `helper_functions.R` and are documented with roxygen style annotations.

## Tech stack

R and the Bioconductor ecosystem: GEOquery, SummarizedExperiment, edgeR, EDASeq, clusterProfiler, enrichplot, org.Hs.eg.db and gprofiler2, together with ggplot2, pheatmap, ggraph and ggkegg for visualization.
